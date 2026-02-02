import xarray as xr
import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

# 1. Load data
ds = xr.open_dataset("all_sites_velocity_subcubes.nc")
site_df = pd.read_csv("points_across_glacier.csv", encoding="utf-8-sig")

# 2. Loop over sites in the xarray
for i, site in enumerate(ds["site"].values):
    # Get initial x, y for this site (assuming 1D coordinates)
    # If x and y are coordinates in ds, you might have ds['x'][i], ds['y'][i]
    # If not, use the CSV to look up the coordinates
    site_name = str(site)
    row = site_df[site_df["Site"] == site_name]
    if row.empty:
        print(f"Site {site_name} not found in CSV, skipping.")
        continue
    x0 = float(row["Lon"].iloc[0])
    y0 = float(row["Lat"].iloc[0])

    target_date = np.datetime64("2023-05-13")
    mid_dates = ds["mid_date"].values
    date_idx = np.abs(mid_dates - target_date).argmin()

    df = ds.sel(site=site)
    dm = df.isel(mid_date=date_idx)

    x_grid = dm["x"].values
    y_grid = dm["y"].values

    # Create interpolators for vx and vy at that time
    vx_interp = RegularGridInterpolator(
        (x_grid, y_grid), dm["vx"].values.T, bounds_error=False, fill_value=np.nan
    )
    vy_interp = RegularGridInterpolator(
        (x_grid, y_grid), dm["vy"].values.T, bounds_error=False, fill_value=np.nan
    )

    print(vx_interp([[x0, y0]]))  # Should give interpolated vx at x0, y0

    def velocity_field(t, Z):
        point = [Z[0], Z[1]]  # x, y
        vx = vx_interp(point)
        vy = vy_interp(point)
        return [vx, vy]

    print(f"Calculating flow line for site {site_name} at coordinates ({x0}, {y0})")

    initial_position = [x0, y0]
    t_span = (0, 10)  # Adjust as needed
    solution = solve_ivp(velocity_field, t_span, initial_position, dense_output=True)

    # Example: print or plot the flow line
    print(f"Site {site_name} flow line X:", solution.y[0])
    print(f"Site {site_name} flow line Y:", solution.y[1])
    # You can plot or save the results as needed
