import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer

ds = xr.open_dataset("resampled_velocity_timeseries.nc")
site_df = pd.read_csv("May2025_fieldsites.csv", encoding='utf-8-sig')

to_meters = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
to_latlon = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)

start_lons = []
start_lats = []
end_lons = []
end_lats = []
lons = []
lats = []
final_sites = []


for idx, site in enumerate(ds['site'].values):
    site_name = str(site)
    df_site = ds.sel(site=site_name)

    row = site_df[site_df['Site'] == site_name]
    if row.empty:
        print(f"Site {site_name} not found in CSV, skipping.")
        continue

    lon0 = float(row['Lon'].iloc[0])
    lat0 = float(row['Lat'].iloc[0])
    x0, y0 = to_meters.transform(lon0, lat0)

    times = pd.to_datetime(df_site['mid_date'].values)
    start_date = pd.Timestamp("2024-05-13")
    end_date = pd.Timestamp("2024-09-04")
    mask = (times >= start_date) & (times <= end_date)
    times = times[mask]
    df_site = df_site.sel(mid_date=times)

    vx = df_site['vx'].values
    vy = df_site['vy'].values

    dt_days = 12
    vx_step = vx * (dt_days / 365.25)
    vy_step = vy * (dt_days / 365.25)

    x = [x0]
    y = [y0]

    for t in range(1, len(times)):
        x_new = x[-1] + vx_step[t - 1]
        y_new = y[-1] + vy_step[t - 1]
        x.append(x_new)
        y.append(y_new)

    x_end = x[-1]
    y_end = y[-1]

    lon_start, lat_start = to_latlon.transform(x0, y0)
    lon_all, lat_all = to_latlon.transform(x, y)
    lon_final, lat_final = to_latlon.transform(x_end, y_end)

    start_lons.append(lon_start)
    start_lats.append(lat_start)
    lons.append(lon_all)
    lats.append(lat_all)
    end_lons.append(lon_final)
    end_lats.append(lat_final)
    final_sites.append(site_name)

# Save final projected positions to CSV
final_df = pd.DataFrame({
    'Site': ds['site'].values,
    'Lat': end_lats,
    'Lon': end_lons
})
final_df.to_csv("Projected_positions.csv", index=False)


plt.figure(figsize=(10, 6))
plt.scatter(start_lons, start_lats, color='blue', label='Initial Position')
plt.scatter(end_lons, end_lats, color='red', label='Final Position')
plt.scatter(lons, lats , color='gray', alpha=0.5, label='Projected Path', s=1.5)

for lon, lat, name in zip(start_lons, start_lats, final_sites):
    plt.text(lon, lat, name, fontsize=12, ha='left', va='bottom')


plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.legend()
plt.title("Initial and projected site positions using ITS_LIVE Velocity Data (May 13, 2024 - Sep 4, 2024)") 

plt.grid()
plt.tight_layout()

plt.savefig("Projected_Flow_Paths.png", dpi=300)

plt.show()