import multiprocessing
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import scipy.sparse
import tqdm
import xarray

import datacube_tools


def single_point_inversion(
    xy, lt, sat_filt=None, start_date=None, stop_date=None, baseline_mask=(None, None), return_data=False,
):
    """
    Invert for a velocity timeseries at a single point. Regularize in space.
    Point must be supplied in EPSG 3413 coordinates.
    Accepts restrictions on satellites used (options = ["1A", "1B", "2A", "2B", "4", "5", "7", "8", "9"])
    example: sat_filt = ["1A", "1B"]
    Also accepts restrictions on date range. Must be supplied as ISO 8601 compliant string
    example: start_date = "2018-01-01", stop_date = "2024-01-01"

    Parameters
    ----------
    xy - coordinates of point in EPSG 3413 (x, y)
    lt - time regularization strength
    sat_filt - list of satellites to use (optional)
    start_date - start date for inversion (optional)
    stop_date - stop date for inversion (optional)
    baseline_mask - (int, int) low and high bound for baselines in integer days.
        None means no bound. Example for no low bound and upper bound of 100 days : (None, 100) (optional)
    return_data - option to return itslive data (optional)

    Returns
    ------
    xres - datacube of inverted velocities
    (if return_data)
    xpnt - itslive datacube
    mask - mask baded on sat_filt, start, stop
    """
    # Download point date
    dc = datacube_tools.DATACUBETOOLS()
    xfull, xpnt, xy = dc.get_timeseries_at_point(
        xy,
        "3413",
        variables=[
            "vx",
            "vy",
            "satellite_img1",
            "satellite_img2",
            "acquisition_date_img1",
            "acquisition_date_img2",
        ],
    )
    xpnt = xpnt.sortby("mid_date")

    # Build mask
    mask = build_mask(xpnt, sat_filt, start_date, stop_date, baseline_mask)

    # Apply mask to velocity grid
    vx_mask = xpnt["vx"][mask].to_numpy()[:, np.newaxis, np.newaxis]
    vy_mask = xpnt["vy"][mask].to_numpy()[:, np.newaxis, np.newaxis]

    # Get masked dates, converting to decimal year
    ns_in_year = 1e9 * 60 * 60 * 24 * 365
    epoch_year = 2000
    epoch = np.datetime64(str(epoch_year))
    aq1 = epoch_year + (
        (xpnt["acquisition_date_img1"][mask] - epoch).to_numpy().astype(np.float64)
        / ns_in_year
    )
    aq2 = epoch_year + (
        (xpnt["acquisition_date_img2"][mask] - epoch).to_numpy().astype(np.float64)
        / ns_in_year
    )

    # Get unique dates
    uniq_aq = np.unique(np.append(aq1, aq2))

    # Solution dates
    solu_dates = (uniq_aq[:-1] + uniq_aq[1:]) / 2

    _, result = solve_point(((0, 0), vx_mask, vy_mask, aq1, aq2, solu_dates, lt))

    # Decimal year to datetime64
    s_in_year = 60 * 60 * 24 * 365
    solu_datetime = [None] * len(solu_dates)
    for i in range(len(solu_dates)):
        whole_year = int(solu_dates[i])
        frac_year = solu_dates[i] - whole_year
        solu_datetime[i] = np.datetime64(str(whole_year))
        solu_datetime[i] += np.timedelta64(int(frac_year * s_in_year), "s")

    # Package into xarray
    xres = xarray.Dataset(
        data_vars=dict(
            vx=(["time"], result["vx"]),
            vy=(["time"], result["vy"]),
        ),
        coords=dict(
            x=xpnt.x,
            y=xpnt.y,
            time=solu_datetime,
        ),
        attrs=dict(
            projection=xpnt.projection,
            GDAL_AREA_OR_POINT=xpnt.GDAL_AREA_OR_POINT,
        ),
    )

    if return_data:
        return xres, xpnt, mask
    else:
        return xres


def solve_point(args):
    """
    Inverts for the velocity time series of a single cell.

    Parameters
    ----------
    args - list of arguments for this function, must be this way to use with imap (sorry).
           ((i,j), vx_mask, vy_mask, aq1, aq2, solu_dates, lt)

    Returns
    -------
    (i,j) - grid cell indices
    cell - vx and vy solutions for the grid cell
    """
    # Unpack arguments
    i, j = args[0]
    vx_mask = args[1]
    vy_mask = args[2]
    aq1 = args[3]
    aq2 = args[4]
    solu_dates = args[5]
    lt = args[6]

    # Invert each component
    cell = {}
    for comp, vgrid in [("vx", vx_mask), ("vy", vy_mask)]:
        # Get unique dates at point
        sub_mask = np.logical_not(np.isnan(vgrid[:, i, j]))

        sub_aq1 = aq1[sub_mask]
        sub_aq2 = aq2[sub_mask]
        sub_uniq_aq = np.unique(np.append(sub_aq1, sub_aq2))
        sub_solu_dates = (sub_uniq_aq[:-1] + sub_uniq_aq[1:]) / 2

        # Quit if there are not enough dates for time regularization to work
        if len(sub_solu_dates) < 3:
            cell = {}
            cell["vx"] = np.zeros(len(solu_dates))
            cell["vy"] = np.zeros(len(solu_dates))
            return (args[0], cell)

        # Set up linear system
        A, b = gen_vel_sys(vgrid[sub_mask, i, j], sub_aq1, sub_aq2, sub_uniq_aq)

        # Add time regularization
        treg = gen_time_regu_second_cent(sub_uniq_aq, lt)
        Aaug = np.vstack((A, treg))
        baug = np.append(b, np.zeros(treg.shape[0]))

        # Solve system
        res = np.linalg.lstsq(Aaug, baug, rcond=None)[0]

        # Interpolate to grid dates
        spl = scipy.interpolate.CubicSpline(sub_solu_dates, res)
        cell[comp] = spl(solu_dates)

    return (args[0], cell)


def solve_stencil(args):
    """
    Inverts a 5 point stencil of velocities. Returns the velocity time series of the middle cell.

    Parameters
    ----------
    args - list of arguments for this function, must be this way to use with imap (sorry).
           ((i,j), vx_mask, vy_mask, aq1, aq2, solu_dates, lt, lx, dx)

    Returns
    -------
    (i,j) - grid cell indices
    cell - vx and vy solutions for the grid cell
    """
    # Unpack arguments
    i, j = args[0]
    vx_mask = args[1]
    vy_mask = args[2]
    aq1 = args[3]
    aq2 = args[4]
    solu_dates = args[5]
    lt = args[6]
    lx = args[7]
    dx = args[8]

    cell = {}
    for comp, vgrid in [("vx", vx_mask), ("vy", vy_mask)]:
        # Make 5 point stencil
        pixels = [
            vgrid[:, i, j],
            vgrid[:, i, j - 1],
            vgrid[:, i, j + 1],
            vgrid[:, i - 1, j],
            vgrid[:, i + 1, j],
        ]

        # Get unique dates in 5 point stencil
        not_nans = [np.logical_not(np.isnan(pix)) for pix in pixels]
        sub_mask = np.logical_or.reduce(not_nans)

        sub_aq1 = aq1[sub_mask]
        sub_aq2 = aq2[sub_mask]
        sub_uniq_aq = np.unique(np.append(sub_aq1, sub_aq2))
        sub_solu_dates = (sub_uniq_aq[:-1] + sub_uniq_aq[1:]) / 2

        # Quit if there are not enough dates for time regularization to work
        if len(sub_solu_dates) < 3:
            cell = {}
            cell["vx"] = np.zeros(len(solu_dates))
            cell["vx"] = np.zeros(len(solu_dates))
            return (args[0], cell)

        # Set up linear systems
        As = []
        bs = [None] * len(pixels)
        for k, v in enumerate(pixels):
            px_mask = np.logical_not(np.isnan(v))
            chunk = [None] * len(pixels)
            chunk[k], bs[k] = gen_vel_sys(
                v[px_mask], aq1[px_mask], aq2[px_mask], sub_uniq_aq
            )
            As.append(chunk)

        # Time regularization
        tregs = []
        for k in range(len(pixels)):
            chunk = [None] * len(pixels)
            chunk[k] = gen_time_regu_second_cent(sub_uniq_aq, lt)
            tregs.append(chunk)

        # Space regularization
        xreg = gen_space_regu_laplacian(sub_uniq_aq, dx, lx)

        # Assemble system
        bs = np.vstack(bs)
        A_sparse = scipy.sparse.bmat(As)
        treg_sparse = scipy.sparse.bmat(tregs)
        Aaug_sparse = scipy.sparse.vstack((A_sparse, treg_sparse, xreg))
        baug = np.vstack(
            (
                bs,
                np.zeros((treg_sparse.shape[0], 1)),
                np.zeros((xreg.shape[0], 1)),
            )
        )[:, 0]

        # Solve system
        res = scipy.sparse.linalg.lsqr(Aaug_sparse, baug)[0][: len(sub_solu_dates)]

        # Interpolate to grid dates
        spl = scipy.interpolate.CubicSpline(sub_solu_dates, res)
        cell[comp] = spl(solu_dates)

    return (args[0], cell)


def grid_inversion(
    xy,
    half_dist,
    lt,
    lx,
    sat_filt=None,
    start_date=None,
    stop_date=None,
    baseline_mask=(None, None),
    return_data=False,
    pbar=False,
    ncpu=1,
):
    """
    Invert for a velocity timeseries for each cell in a velocity grid. Regularize in time and space.
    Use multiprocessing to speed up the computation.
    Center point of grid must be supplied in EPSG 3413 coordinates.
    Accepts restrictions on satellites used (options = ["1A", "1B", "2A", "2B", "4", "5", "7", "8", "9"])
    example: sat_filt = ["1A", "1B"]
    Also accepts restrictions on date range. Must be supplied as ISO 8601 compliant string
    example: start_date = "2018-01-01", stop_date = "2024-01-01"

    Parameters
    ----------
    xy - coordinates of grid center in EPSG 3413 (x, y)
    half_dist - half width of grid in meters
    lt - time regularization strength
    lx - space regularization strength
    sat_filt - list of satellites to use (optional)
    start_date - start date for inversion (optional)
    stop_date - stop date for inversion (optional)
    baseline_mask - (int, int) low and high bound for baselines in integer days.
        None means no bound. Example for no low bound and upper bound of 100 days : (None, 100) (optional)
    return_data - option to return itslive data
    pbar - option to show progress bar
    ncpu - number of processes to use (optional, default 1)

    Returns
    ------
    xres - inverted velocity datacube
    (if return_data)
    xsub - itslive datacube
    mask - mask baded on sat_filt, start, stop
    """
    # Download point date
    dc = datacube_tools.DATACUBETOOLS()
    xfull, xsub, xy = dc.get_subcube_around_point(
        xy,
        "3413",
        half_distance=half_dist,
        variables=[
            "vx",
            "vy",
            "satellite_img1",
            "satellite_img2",
            "acquisition_date_img1",
            "acquisition_date_img2",
        ],
    )
    xsub = xsub.sortby("mid_date")

    if lx != 0:
        if np.any(np.array(xsub["vx"].shape[1:]) < 3):
            print("Grid not large enough for 5 point stencil spatial regularization")
            exit()

    # Build mask
    mask = build_mask(xsub, sat_filt, start_date, stop_date, baseline_mask)

    # Get masked dates, converting to decimal year
    ns_in_year = 1e9 * 60 * 60 * 24 * 365
    epoch_year = 2015
    epoch = np.datetime64(str(epoch_year))
    aq1 = epoch_year + (
        (xsub["acquisition_date_img1"][mask] - epoch).to_numpy().astype(np.float64)
        / ns_in_year
    )
    aq2 = epoch_year + (
        (xsub["acquisition_date_img2"][mask] - epoch).to_numpy().astype(np.float64)
        / ns_in_year
    )

    # Get unique dates
    uniq_aq = np.unique(np.append(aq1, aq2))

    # Solution dates
    solu_dates = (uniq_aq[:-1] + uniq_aq[1:]) / 2

    # Apply mask to velocity grid
    vx_mask = xsub["vx"][mask, :, :].to_numpy()
    vy_mask = xsub["vy"][mask, :, :].to_numpy()

    # Invert each component
    result = {
        "vx": np.zeros((len(solu_dates), vx_mask.shape[1], vx_mask.shape[2])),
        "vy": np.zeros((len(solu_dates), vy_mask.shape[1], vy_mask.shape[2])),
    }

    # Make list of pixels
    argslist = []
    if lx == 0:
        dx = 0
        for i in range(1, vx_mask.shape[1] - 1):
            for j in range(1, vx_mask.shape[2] - 1):
                argslist.append(((i, j), vx_mask, vy_mask, aq1, aq2, solu_dates, lt))
    else:
        dx = np.diff(xsub.x)[0]
        for i in range(1, vx_mask.shape[1] - 1):
            for j in range(1, vx_mask.shape[2] - 1):
                argslist.append(
                    ((i, j), vx_mask, vy_mask, aq1, aq2, solu_dates, lt, lx, dx)
                )

    # Run inversions
    with multiprocessing.Pool(ncpu) as p:
        if lx == 0:
            cells = list(
                tqdm.tqdm(
                    p.imap(solve_point, argslist), total=len(argslist), disable=not pbar
                )
            )
        else:
            cells = list(
                tqdm.tqdm(
                    p.imap(solve_stencil, argslist),
                    total=len(argslist),
                    disable=not pbar,
                )
            )

    # Unpack results
    for cell in cells:
        (i, j), res = cell
        for comp, v in res.items():
            result[comp][:, i, j] = v[:]

    # Decimal year to datetime64
    s_in_year = 60 * 60 * 24 * 365
    solu_datetime = [None] * len(solu_dates)
    for i in range(len(solu_dates)):
        whole_year = int(solu_dates[i])
        frac_year = solu_dates[i] - whole_year
        solu_datetime[i] = np.datetime64(str(whole_year))
        solu_datetime[i] += np.timedelta64(int(frac_year * s_in_year), "s")

    # Package into xarray
    xres = xarray.Dataset(
        data_vars=dict(
            vx=(["time", "y", "x"], result["vx"]),
            vy=(["time", "y", "x"], result["vy"]),
        ),
        coords=dict(
            x=xsub.x,
            y=xsub.y,
            time=solu_datetime,
        ),
        attrs=dict(
            projection=xsub.projection,
            GDAL_AREA_OR_POINT=xsub.GDAL_AREA_OR_POINT,
        ),
    )

    if return_data:
        return xres, xsub, mask
    else:
        return xres


def build_mask(dcube, sat_filt, start_date, stop_date, baseline_mask):
    """
    Build velocity mask for inversion based on NaNs, satellite choices, and a date range

    Parameters
    ----------
    dcube - velocity datecube
    sat_filt - list of satellites to use
    start_date - start date for inversion
    stop_date - stop date for inversion
    baseline_mask - tuple specifying high and low baseline bounds in integer days

    Returns
    ------
    mask - Boolean mask for time axis of velocity datacube
    """
    # Remove all-nan timesteps, use vx
    if len(dcube.vx.shape) == 3:
        # supress stupid nanmean warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mask = np.logical_not(
                np.isnan(np.nanmean(np.nanmean(dcube.vx, axis=1), axis=1))
            )
    elif len(dcube.vx.shape) == 1:
        mask = np.logical_not(np.isnan(dcube.vx))

    # Satellite filter
    if sat_filt is not None:
        mask = np.logical_and(
            mask,
            np.logical_and(
                np.in1d(dcube.satellite_img1, sat_filt),
                np.in1d(dcube.satellite_img2, sat_filt),
            ),
        )

    # Start date filter
    if start_date is not None:
        tstart = np.datetime64(start_date)
        mask = np.logical_and(mask, dcube.acquisition_date_img1 >= tstart)

    # Stop date filter
    if stop_date is not None:
        tstop = np.datetime64(stop_date)
        mask = np.logical_and(mask, dcube.acquisition_date_img2 <= tstop)

    # Baseline filter
    baseline = (dcube.acquisition_date_img2 - dcube.acquisition_date_img1).dt.days
    
    if(baseline_mask[0] is not None):
        mask = np.logical_and(mask, baseline > baseline_mask[0])

    if(baseline_mask[1] is not None):
        mask = np.logical_and(mask, baseline < baseline_mask[1])
    
    return mask


def gen_vel_sys(vel, aq1, aq2, uniq_aq):
    """
    Generate matrix representation of system of velocity equations

    Parameters
    ----------
    vel - velocity timeseries (m per year)
    aq1 - First aquisition date for each velocity in timeseries (decimal year)
    aq2 - Second acquisiton date for each velocity in timeseries (decimal year)
    uniq_aq - Time of unique aquisitions (decimal year)

    Returns
    -------
    A, b - Matrices representing system of velocity equations
    """
    uniq_aq = np.unique(uniq_aq)
    diff_aq = np.diff(uniq_aq)

    A = np.zeros((len(vel), len(diff_aq)), dtype=np.float64)
    b = np.zeros((len(vel), 1), dtype=np.float64)

    for i, v in enumerate(vel):
        j0 = np.where(aq1[i] == uniq_aq)[0][0]
        j1 = np.where(aq2[i] == uniq_aq)[0][0]
        A[i, j0:j1] = diff_aq[j0:j1]
        b[i] = vel[i] * (aq2[i] - aq1[i])

    return A, b


def gen_time_regu_first_fwd(uniq_aq, l):
    """
    Generate time regularization matrix - first order forward difference.
    Handles variable time steps.

    Parameters
    ----------
    uniq_aq - unique acquisition dates (decimal year)
    l - regularization strength

    Returns
    -------
    regt - regularization matrix
    """
    uniq_aq = np.sort(uniq_aq)
    diff_aq = np.diff(uniq_aq)
    treg = np.zeros((len(diff_aq) - 1, len(diff_aq)))
    for i in range(len(diff_aq) - 1):
        treg[i, i] = l / diff_aq[i]
        treg[i, i + 1] = -l / diff_aq[i]

    return regt


def gen_time_regu_second_cent(uniq_aq, l):
    """
    Generate time regularization matrix - second order central difference.
    Handles variable timesteps.

    Parameters
    ----------
    uniq_aq - unique acquisition dates (decimal year)
    l - regularization strength

    Returns
    -------
    regt - regularization matrix
    """
    uniq_aq = np.sort(uniq_aq)
    diff_aq = np.diff(uniq_aq)
    treg = np.zeros((len(diff_aq) - 2, len(diff_aq)))
    for i in range(len(diff_aq) - 2):
        treg[i, i] = l * 2 / (diff_aq[i] * (diff_aq[i] + diff_aq[i + 1]))
        treg[i, i + 1] = -l * 2 / (diff_aq[i] * diff_aq[i + 1])
        treg[i, i + 2] = l * 2 / (diff_aq[i + 1] * (diff_aq[i] + diff_aq[i + 1]))

    return treg


def gen_space_regu_laplacian(uniq_aq, dx, l):
    """
    Generate space regularization matrix - five point stencil laplacian.
                                                       4
    Assumes cells are in this order in the matrix:   2 1 3
                                                       5
    Parameters
    ----------
    uniq_aq - unique acquisition dates (decimal year)
    dx - spatial discretization (m)
    l - regularization strength

    Returns
    -------
    regt - regularization matrix
    """
    uniq_aq = np.sort(uniq_aq)
    diff_aq = np.diff(uniq_aq)
    xreg = np.zeros((len(diff_aq), len(diff_aq) * 5))
    for i in range(len(diff_aq)):
        xreg[i, i] = -4 * l / (dx**2)
        xreg[i, i + len(diff_aq)] = l / (dx**2)
        xreg[i, i + 2 * len(diff_aq)] = l / (dx**2)
        xreg[i, i + 3 * len(diff_aq)] = l / (dx**2)
        xreg[i, i + 4 * len(diff_aq)] = l / (dx**2)

    return xreg
