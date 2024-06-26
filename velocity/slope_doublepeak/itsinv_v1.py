import numpy as np
import datacube_tools
import scipy.sparse
import tqdm
import warnings


def single_point_inversion(
    xy, lt, sat_filt=None, start_date=None, stop_date=None, return_data=False
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
    return_data - option to return itslive data

    Returns
    ------
    solu_dates - mid dates of average velocity solutions
    result - dictionary of average velocity solutions
    (if return_data)
    mid_dates - mid dates of itslive data
    itsl - dictionary of itslive data
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
    mask = build_mask(xpnt, sat_filt, start_date, stop_date)

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

    # Invert each component
    result = {}
    for comp in ["vx", "vy"]:
        # Set up linear system
        A, b = gen_vel_sys(xpnt[comp][mask].to_numpy(), aq1, aq2, uniq_aq)

        # Add time regularization
        treg = gen_time_regu_second_cent(uniq_aq, lt)
        Aaug = np.vstack((A, treg))
        baug = np.append(b, np.zeros(treg.shape[0]))

        # Solve system
        result[comp] = np.linalg.lstsq(Aaug, baug, rcond=None)[0]

    if return_data:
        mid_dates = (aq1 + aq2) / 2
        itsl = {}
        for comp in ["vx", "vy"]:
            itsl[comp] = xpnt[comp][mask].to_numpy()
        return solu_dates, result, mid_dates, itsl
    else:
        return solu_dates, result


def grid_inversion(
    xy,
    half_dist,
    lt,
    lx,
    sat_filt=None,
    start_date=None,
    stop_date=None,
    return_data=False,
    pbar=False,
):
    """
    Invert for a velocity timeseries for each cell in a velocity grid. Regularize in time and space.
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
    return_data - option to return itslive data
    pbar - option to show progress bar

    Returns
    ------
    solu_dates - mid dates of average velocity solutions
    result - dictionary of average velocity solutions
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

    # Build mask
    mask = build_mask(xsub, sat_filt, start_date, stop_date)

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
    # Progress bar
    with tqdm.tqdm(
        total=(vx_mask.shape[1] - 2) * (vx_mask.shape[2] - 2) * 2, disable=not pbar
    ) as pbar:
        # Loops over grid
        for i in range(1, vx_mask.shape[1] - 1):
            for j in range(1, vx_mask.shape[2] - 1):
                for comp, vgrid in [("vx", vx_mask), ("vy", vy_mask)]:
                    # Set up linear systems
                    As = []
                    bs = [None] * 5
                    for k, v in enumerate(
                        [
                            vgrid[:, i, j],
                            vgrid[:, i, j - 1],
                            vgrid[:, i, j + 1],
                            vgrid[:, i - 1, j],
                            vgrid[:, i + 1, j],
                        ]
                    ):
                        px_mask = np.logical_not(np.isnan(v))
                        chunk = [None] * 5
                        chunk[k], bs[k] = gen_vel_sys(
                            v[px_mask], aq1[px_mask], aq2[px_mask], uniq_aq
                        )
                        As.append(chunk)

                    # Time regularization
                    tregs = []
                    for k in range(5):
                        chunk = [None] * 5
                        chunk[k] = gen_time_regu_second_cent(uniq_aq, lt)
                        tregs.append(chunk)

                    # Space regularization
                    dx = np.diff(xsub.x)[0]
                    xreg = gen_space_regu_laplacian(uniq_aq, dx, lx)

                    # Assemble system
                    A_sparse = scipy.sparse.bmat(As)
                    treg_sparse = scipy.sparse.bmat(tregs)
                    Aaug_sparse = scipy.sparse.vstack((A_sparse, treg_sparse, xreg))

                    bs = np.vstack(bs)
                    baug = np.vstack(
                        (
                            bs,
                            np.zeros((treg_sparse.shape[0], 1)),
                            np.zeros((xreg.shape[0], 1)),
                        )
                    )[:, 0]

                    # Solve system
                    result[comp][:, i, j] = scipy.sparse.linalg.lsqr(Aaug_sparse, baug)[
                        0
                    ][: len(solu_dates)]

                    pbar.update(1)

    if return_data:
        return solu_dates, result, xsub, mask
    else:
        return solu_dates, result


def build_mask(dcube, sat_filt, start_date, stop_date):
    """
    Build velocity mask for inversion based on NaNs, satellite choices, and a date range

    Parameters
    ----------
    dcube - velocity datecube
    sat_filt - list of satellites to use
    start_date - start date for inversion
    stop_date - stop date for inversion

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
        xreg[i, i] = -4 * l / (dx ** 2)
        xreg[i, i + len(diff_aq)] = l / (dx ** 2)
        xreg[i, i + 2 * len(diff_aq)] = l / (dx ** 2)
        xreg[i, i + 3 * len(diff_aq)] = l / (dx ** 2)
        xreg[i, i + 4 * len(diff_aq)] = l / (dx ** 2)

    return xreg
