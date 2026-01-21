import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm

##### Load the velocity time series dataset #####
df = xr.open_dataset("velocity_timeseries.nc")
print("Dataset loaded successfully.")

##### Filter out long baselines #####
baseline_threshold = pd.Timedelta(days=20)
df = df.where(df['date_dt'] < baseline_threshold, drop=True)
df = df.dropna(dim='mid_date', how='any')
df = df.sortby('mid_date')

##### Drop unnecessary data from time series ####
start_date = np.datetime64("2014-10-01")
ds = df.where(df['mid_date'] > start_date, drop=True)

###### LOWESS smoothing for v, vx, and vy ######
dates = pd.to_datetime(ds['mid_date'].values)
v_lowess = np.full(ds['v'].shape, np.nan)
vx_lowess = np.full(ds['vx'].shape, np.nan)
vy_lowess = np.full(ds['vy'].shape, np.nan)

for i, site in enumerate(ds['site'].values):
    # Smooth v
    v = ds['v'][i].values
    mask_v = ~np.isnan(v)
    if np.sum(mask_v) >= 5:
        t_index = np.array([(d - dates[mask_v][0]).days for d in dates[mask_v]])
        v_smooth = sm.nonparametric.lowess(v[mask_v], t_index, frac=1/25, return_sorted=False)
        v_lowess[i, mask_v] = v_smooth

    # Smooth vx
    vx = ds['vx'][i].values
    mask_vx = ~np.isnan(vx)
    if np.sum(mask_vx) >= 5:
        t_index_vx = np.array([(d - dates[mask_vx][0]).days for d in dates[mask_vx]])
        vx_smooth = sm.nonparametric.lowess(vx[mask_vx], t_index_vx, frac=1/25, return_sorted=False)
        vx_lowess[i, mask_vx] = vx_smooth

    # Smooth vy
    vy = ds['vy'][i].values
    mask_vy = ~np.isnan(vy)
    if np.sum(mask_vy) >= 5:
        t_index_vy = np.array([(d - dates[mask_vy][0]).days for d in dates[mask_vy]])
        vy_smooth = sm.nonparametric.lowess(vy[mask_vy], t_index_vy, frac=1/25, return_sorted=False)
        vy_lowess[i, mask_vy] = vy_smooth

# Add smoothed data to dataset
ds['vs'] = (('site', 'mid_date'), v_lowess)
ds['vxs'] = (('site', 'mid_date'), vx_lowess)
ds['vys'] = (('site', 'mid_date'), vy_lowess)

###### Resample and interpolate ######
ds_resample = ds.resample(mid_date='12D').mean()
ds_resample = ds_resample[['vs', 'vxs', 'vys']]
ds_resample = ds_resample.interpolate_na(dim='mid_date', method='linear')
ds_resample = ds_resample.ffill('mid_date').bfill('mid_date')

# Rename and save to final dataset
dr = xr.Dataset(
    {
        'v': ds_resample['vs'],
        'vx': ds_resample['vxs'],
        'vy': ds_resample['vys']
    },
    coords={
        'mid_date': ds_resample['mid_date'],
        'site': ds.site
    }
)
dr.to_netcdf("resampled_velocity_timeseries.nc")
print("Smoothed and resampled dataset saved.")

###### Plotting examples ######

# vx plot
plt.figure(figsize=(12, 4))
for i in range(ds.dims['site']):
    site = ds.site[i].values
    plt.plot(dr['mid_date'], dr['vx'][:, i], label=f'{site} vx')
plt.xlabel('Date')
plt.ylabel('vx (m/yr)')
plt.title('Smoothed & Resampled vx')
plt.legend()
plt.tight_layout()
plt.show()

# vy plot
plt.figure(figsize=(12, 4))
for i in range(ds.dims['site']):
    site = ds.site[i].values
    plt.plot(dr['mid_date'], dr['vy'][:, i], label=f'{site} vy')
plt.xlabel('Date')
plt.ylabel('vy (m/yr)')
plt.title('Smoothed & Resampled vy')
plt.legend()
plt.tight_layout()
plt.show()

# v (magnitude) plot
plt.figure(figsize=(12, 4))
for i in range(ds.dims['site']):
    site = ds.site[i].values
    plt.plot(dr['mid_date'], dr['v'][:, i], label=f'{site}')
plt.xlabel('Date')
plt.ylabel('v (m/yr)')
plt.title('Smoothed and resampled speed (v)')
plt.legend()
plt.tight_layout()
plt.savefig("smoothed_resampled_v_magnitude.png")
plt.show()

# Compute v from vx and vy
v_calc = np.sqrt(dr['vx']**2 + dr['vy']**2)

# Plot comparison of original v vs computed v from vx and vy
plt.figure(figsize=(12, 4))
for i in range(dr.dims['site']):
    site = dr.site[i].values
    plt.plot(dr['mid_date'], dr['v'][:, i], label=f'{site} v (original)', linestyle='-')
    plt.plot(dr['mid_date'], v_calc[:, i], label=f'{site} v (from vx, vy)', linestyle='--')

plt.xlabel('Date')
plt.ylabel('v (m/yr)')
plt.title('Comparison: Original vs Computed v from vx, vy')
plt.legend()
plt.tight_layout()
plt.show()
