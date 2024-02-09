import argparse
import os
import warnings

import cartopy
import matplotlib.pyplot as plt
import numpy as np
import scipy
import tqdm
import xarray


def cli():
    parser = argparse.ArgumentParser(
        prog="make_movie_frames.py", description="Make frames for a velocity movie"
    )
    parser.add_argument(
        "datacube", help="Path to inversion result xarray (saved as netcdf)"
    )
    return parser.parse_args()


def main():
    args = cli()

    vel = xarray.open_dataset(args.datacube)

    # Interpolate to daily velocities
    s_in_day = 60 * 60 * 24
    epoch = 2015
    epoch = np.datetime64(str(epoch), "s")

    # datetime64 bounds and steps
    t0 = np.min(vel["time"].to_numpy())
    t1 = np.max(vel["time"].to_numpy())
    dt = np.timedelta64(1, "D")
    nstep = int((t1 - t0) / dt)
    tsteps = np.array([t0 + dt * i for i in range(nstep)])

    # decimal day bounds and steps
    # t0d = ((t0.astype('datetime64[s]')-epoch).astype(np.float32)/s_in_day)
    # t1d = ((t1.astype('datetime64[s]')-epoch).astype(np.float32)/s_in_day)
    # dtd = 1
    # tstepsd = np.array([t0d+dtd*i for i in range(nstep)])

    # decimal days of data points
    # tvd = ((vel["time"].to_numpy().astype('datetime64[s]')-epoch).astype(np.float32)/s_in_day)

    # Make interpolated dataset
    # vxi = np.zeros((nstep, vel["vx"].shape[1], vel["vx"].shape[2]))
    # vyi = np.zeros((nstep, vel["vy"].shape[1], vel["vy"].shape[2]))

    # for i in tqdm.tqdm(range(vel["vx"].shape[1]), desc="Interpolating to daily velocities"):
    #    for j in range(vel["vx"].shape[2]):
    #        splinex = scipy.interpolate.UnivariateSpline(tvd, vel["vx"][:, i, j], s=len(tvd)//2)
    #        spliney = scipy.interpolate.UnivariateSpline(tvd, vel["vy"][:, i, j], s=len(tvd)//2)
    #        vxi[:, i, j] = splinex(tstepsd)
    #        vyi[:, i, j] = spliney(tstepsd)

    # Cubic spline interpolation
    splinex = scipy.interpolate.CubicSpline(vel["time"], vel["vx"])
    spliney = scipy.interpolate.CubicSpline(vel["time"], vel["vy"])

    vxi = splinex(tsteps)
    vyi = spliney(tsteps)

    # Generate frames
    prj = cartopy.crs.epsg("3338")
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=prj)
    grid_extent = (vel["x"].min(), vel["x"].max(), vel["y"].min(), vel["y"].max())

    os.system("mkdir -p frames")

    for i, t in tqdm.tqdm(enumerate(tsteps), total=len(tsteps)):
        ax.set_extent([-139.55, -139.31, 59.98, 60.1], crs=cartopy.crs.PlateCarree())
        vx = vxi[i, :, :]
        vy = vyi[i, :, :]
        v = np.sqrt(vx**2 + vy**2)

        # Mask out unreasonable hubbard velocities
        mask = v > 10000
        v[mask] = np.nan
        vx[mask] = 0
        vy[mask] = 0

        date = np.datetime_as_string(t, unit="D")
        ax.set_title(date)
        ax.imshow(
            v,
            transform=cartopy.crs.epsg("3413"),
            cmap="jet",
            extent=grid_extent,
            vmin=0,
            vmax=3500,
        )
        vecds = 3
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ax.quiver(
                vel["x"][::vecds],
                vel["y"][::vecds],
                vx[::vecds, ::vecds],
                vy[::vecds, ::vecds],
                transform=cartopy.crs.epsg("3413"),
                scale=5e4,
            )

        fig.savefig("frames/%s.png" % date, bbox_inches="tight")
        ax.cla()


main()
