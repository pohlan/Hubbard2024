import numpy as np
import datacube_tools
import csv
from velocity_widget import ITSLIVE
velocity_widget = ITSLIVE()
dc = datacube_tools.DATACUBETOOLS()
import xarray as xr
import glob
import os
import shutil

site_data = []
with open("points_across_glacier.csv", newline='',encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['lat'] and row['lon']:  # Only add rows with both values present
            site_data.append(row)

lats = [float(row['lat']) for row in site_data]
lons = [float(row['lon']) for row in site_data]
sites = [row['Site'] for row in site_data]

def reduce_cube(ds):
    # Filter by start date
    start_date = np.datetime64("2023-01-01")
    ds = ds.sel(mid_date=ds['mid_date'] >= start_date)

    # Filter out long baselines
    idx_baselines = np.where(ds['date_dt'].values < np.timedelta64(20, 'D'))[0]
    ds = ds.isel(mid_date=idx_baselines)

    # Drop unnecessary variables
    keep_vars = ['vx', 'vy', 'mid_date']
    drop_vars = [var for var in ds.data_vars if var not in keep_vars]
    ds = ds.drop_vars(drop_vars)

    return ds


for lat, lon, site in zip(lats, lons, sites):
    print(f"Querying and processing {site}...")
    # Query the cube for this site
    cube = dc.get_subcube_around_point(
        point_xy=(lon, lat),
        point_epsg_str=4326,
        half_distance=1500.0,
        variables=["vx", "vy", "mid_date", "date_dt"]
    )
    ds = cube[0]

    # Reduce the cube
    ds_reduced = reduce_cube(ds)

    # Optionally add a site dimension
    ds_reduced = ds_reduced.expand_dims(dim={"site": [site]})

    # Save to NetCDF (or Zarr)
    filename = f"{site}_velocity_subcube.nc"
    print(f"Saving {filename} ...")
    ds_reduced.to_netcdf(filename)

# Make sure the folder exists
os.makedirs("site_subcubes", exist_ok=True)

# Move all matching files to the folder
for f in glob.glob("*_velocity_subcube.nc"):
    shutil.move(f, os.path.join("site_subcubes", os.path.basename(f)))

# Path to your subcube files
folder = "./site_subcubes"

file_list = sorted(glob.glob(os.path.join(folder, "*_velocity_subcube.nc")))

# Load all datasets into a list
datasets = [xr.open_dataset(f) for f in file_list]

# Concatenate along the 'site' dimension
combined = xr.concat(datasets, dim='site')

# Optionally, save the combined dataset
combined.to_netcdf("all_sites_velocity_subcubes.nc")

print("Combined dataset shape:", combined.dims)