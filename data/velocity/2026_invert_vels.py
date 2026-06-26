"""
Grid inversion of ITS_LIVE velocities (converted from invert_annual_velocity.ipynb).

Flow:
  1. Configure the inversion (point, half_dist, sats, dates, regularization).
  2. Plot the cube FOOTPRINT (no download) and wait for confirmation.
  3. Only on confirmation: download the subcube + run the inversion (the heavy part).
  4. Save the result and show diagnostic plots.

Run from this folder so `import itsinv` resolves:
    cd ~/Desktop/PHD/Paper_Hubbard/Hubbard2024-main/velocity
    python invert_annual_velocity.py          # asks before the heavy run
    python invert_annual_velocity.py --yes     # skip the confirmation gate
"""

import argparse

import numpy as np
import pandas as pd
import xarray
import matplotlib.pyplot as plt

# Module names start with a digit, so they can't use plain `import` — load via
# importlib and bind to the short names the rest of this file uses.
import importlib
itsinv = importlib.import_module("2026_itsinv")
datacube_tools = importlib.import_module("2026_datacube_tools")

# =============================================================================
# CONFIG
# =============================================================================
XY = (-3306261, 254285)          # center point, EPSG:3413 (x, y) in metres
NAME = "Test2"
LIST_SATS = ["1A", "1B"]         # satellites to use (or None for all)

HALF_DIST = 10000                # half-width of the cube box, metres
LT = 1e-4                        # time regularization strength
LX = 10                          # space regularization strength

START_DATE = pd.to_datetime("2016-01-01")                # e.g. pd.to_datetime("2016-01-01") or None
STOP_DATE = pd.to_datetime("2025-01-01")                 # e.g. pd.to_datetime("2025-01-01") or None

OUTPUT = f"{NAME}_velocities.nc"
EPSG = "3413"
ITSLIVE_RES_M = 120              # nominal ITS_LIVE pixel size, for the size estimate
PREVIEW_YEAR = 2022              # the single year downloaded for the footprint preview


# =============================================================================
# FOOTPRINT PREVIEW (downloads ONLY 1 year, plots the mean speed, then confirm)
# =============================================================================
def _open_subcube_lazy(xy, half_dist):
    """Open the ITS_LIVE cube lazily and return the spatial box around `xy`
    WITHOUT loading the full time series (unlike grid_inversion's loader)."""
    dc = datacube_tools.DATACUBETOOLS()
    cube_feature, (px, py) = dc.find_datacube_catalog_entry_for_point(xy, EPSG)
    url = (
        cube_feature["properties"]["zarr_url"]
        .replace("https:", "s3:")
        .replace("http:", "s3:")
        .replace(".s3.amazonaws.com", "")
    )
    ds = xarray.open_dataset(url, engine="zarr", storage_options={"anon": True})
    lx, ly = ds.x, ds.y
    sub = ds[["vx", "vy"]].sel(
        x=lx[(lx > px - half_dist) & (lx < px + half_dist)],
        y=ly[(ly > py - half_dist) & (ly < py + half_dist)],
    )
    return sub, (px, py)


def plot_footprint(xy, half_dist, name, year=PREVIEW_YEAR, save="footprint_preview.png"):
    """Download a single year of velocity, plot the mean speed over the cube
    footprint, and show it so the region can be confirmed before the heavy
    download + inversion. Only that year's chunks are read into memory."""
    side_km = 2 * half_dist / 1000
    t0 = f"{year}-01-01"
    t1 = f"{year}-12-31"
    print(f"Preview: downloading {year} velocities for the {name} footprint "
          f"({side_km:.0f} km box)...")

    sub, (px, py) = _open_subcube_lazy(xy, half_dist)
    sub = sub.sortby("mid_date").sel(mid_date=slice(t0, t1))
    n_scenes = sub.sizes.get("mid_date", 0)
    if n_scenes == 0:
        print(f"  No scenes in {year} for this box — try another PREVIEW_YEAR.")
        return
    speed = np.sqrt(sub.vx ** 2 + sub.vy ** 2)
    vmean = speed.mean(dim="mid_date").compute()   # <- only here does data load

    fig, ax = plt.subplots(figsize=(6.5, 6))
    vmean.plot.imshow(
        ax=ax, vmin=0, vmax=2500, cmap="jet",
        cbar_kwargs={"label": "mean speed [m/yr]"},
    )
    ax.plot(px, py, "r*", ms=14, label="center")
    ax.set_aspect("equal")
    ax.set_xlabel(f"EPSG:{EPSG} x  [m]")
    ax.set_ylabel(f"EPSG:{EPSG} y  [m]")
    ax.set_title(f"{name} footprint — mean speed {year}  "
                 f"({n_scenes} scenes, {side_km:.0f} km, half_dist={half_dist} m)")
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(save, dpi=120)
    print(f"  Preview saved to {save}  ({n_scenes} scenes, grid {dict(vmean.sizes)})")
    plt.show()


def confirm(prompt="Proceed with download + inversion? [y/N] "):
    try:
        return input(prompt).strip().lower() in ("y", "yes")
    except EOFError:
        return False


# =============================================================================
# INVERSION
# =============================================================================
def run_inversion():
    print(f"Processing {NAME}  sats={LIST_SATS}  half_dist={HALF_DIST}  "
          f"dates=({START_DATE}, {STOP_DATE})")
    xres, xsub, mask = itsinv.grid_inversion(
        XY,
        half_dist=HALF_DIST,
        lt=LT,
        lx=LX,
        sat_filt=LIST_SATS,
        start_date=START_DATE,
        stop_date=STOP_DATE,
        pbar=True,
        return_data=True,
    )
    xres.to_netcdf(OUTPUT)
    print(f"Saved inversion to {OUTPUT}")
    return xres, xsub, mask


def plot_results(xres, xsub, mask):
    """Diagnostic plots: speed timeseries at the cube center + mean speed map."""
    v = np.sqrt(xres["vx"] ** 2 + xres["vy"] ** 2)
    ny, nx = v.shape[1], v.shape[2]
    cy, cx = ny // 2, nx // 2

    # Center-pixel timeseries: inversion vs. raw data
    md = xsub.acquisition_date_img1 + (
        (xsub.acquisition_date_img2 - xsub.acquisition_date_img1) / 2
    )
    v_data = np.sqrt(xsub.vx ** 2 + xsub.vy ** 2)
    plt.figure()
    plt.plot(md[mask], v_data[mask][:, cy, cx], "k.", label="Data")
    plt.plot(xres.time, v[:, cy, cx], "r-", label="Inversion")
    plt.ylabel("speed [m/yr]")
    plt.title(f"{NAME} center pixel ({cy},{cx})")
    plt.legend()

    # Mean speed map
    plt.figure()
    plt.imshow(np.nanmean(v, axis=0), vmin=0, vmax=2500, cmap="jet")
    plt.colorbar(label="mean speed [m/yr]")
    plt.title(f"{NAME} mean speed")
    plt.show()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-y", "--yes", action="store_true",
                    help="skip the footprint confirmation gate")
    ap.add_argument("--footprint-only", action="store_true",
                    help="only show the footprint preview, then exit")
    args = ap.parse_args()

    plot_footprint(XY, HALF_DIST, NAME)
    if args.footprint_only:
        return
    if not args.yes and not confirm():
        print("Aborted before download/inversion.")
        return

    xres, xsub, mask = run_inversion()
    plot_results(xres, xsub, mask)


if __name__ == "__main__":
    main()
