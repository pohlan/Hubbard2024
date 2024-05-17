import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray
import scipy
from scipy.signal import savgol_filter
import xarray
from tqdm import tqdm
import rasterio as rio
import affine
import os
import elevation
import imageio.v2 as imageio
from matplotlib.colors import BoundaryNorm
import geopandas as gpd
from pyproj import Transformer
import matplotlib.colors as mcolors

########################

# load DEM
    
# Open the GeoTIFF file
with rio.open('ifsar_hubbardDEM_reproj.tif') as src:
    # Read the raster data
    data = src.read(1)  # Assuming a single band image
    height = data.shape[0]
    width = data.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rio.transform.xy(src.transform, rows, cols)
    x_coords= np.array(xs)
    y_coords = np.array(ys)
    
    # Calculate the resolution (pixel size) of the DEM
    resolution = src.res[0]  # Assuming square pixels
    
    # Compute gradients using finite differences
    dz_dx, dz_dy = np.gradient(data, resolution)

    # Calculate slope magnitude
    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    slope = np.degrees(slope_rad)   
    
    # Plot the DEM using matplotlib
    #plt.figure(figsize=(8, 6))
    #plt.imshow(data, cmap='terrain', aspect='auto')
    #plt.colorbar(label= Slope (deg)')
    #plt.title('Digital Elevation Model')
    #plt.grid(True)
    #plt.show()
    
    # Plot the slope using matplotlib
    #plt.figure(figsize=(8, 6))
    #plt.imshow(slope, cmap='terrain', aspect='auto')
    #plt.colorbar(label='Elevation (m)')
    #plt.title('Slope of DEM')
    #plt.grid(True)
    #plt.show()
    
print('finished loading DEM and calculating slope') 

########################
# load centerline

points = pd.read_csv('centerline_points_100m.csv')
x, y =points.X.to_numpy(), points.Y.to_numpy()

#distance from terminus along centerline 
d = np.linspace(100*len(x), 0, len(x))

print('finished loading centerlines')
 
########################
# get slope along centerline

#initialize arrays
slope_centerline = np.zeros((len(x)))
dem_centerline = np.zeros((len(x)))


# find indices along centerline 
for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(x_coords[:, :] - x[i]).argmin(axis = 1)
    target_y_idx= np.abs(y_coords[:, :] - y[i]).argmin(axis = 0)
    
     # Extract the slope value closest to the target coordinates
    slope_centerline[i] = slope[target_y_idx[0], target_x_idx[0]]
    dem_centerline[i] = data[target_y_idx[0], target_x_idx[0]]

# Combine target_x_idx, target_y_idx, and slope into a DataFrame
df1 = pd.DataFrame({'x': x, 'y': y, 'slope': slope_centerline})
df2 = pd.DataFrame({'x': x, 'y': y, 'dem': dem_centerline})
# Write the DataFrame to a CSV file
df1.to_csv('slope_along_centerline.csv')
df2.to_csv('dem_along_centerline.csv')

print('finished calculating slope along centerline')


########################

ds = xarray.load_dataset("../Hubbard_5eminus5.nc")
ds["v"] = np.sqrt(ds.vx**2 + ds.vy**2).fillna(0)
ds["month"] = ds.time.dt.month
ds["year"] = ds.time.dt.year
ds["doy"] = ds.time.dt.dayofyear
xx, yy = np.meshgrid(ds.x, ds.y)

years = [2018, 2019, 2020, 2021]

strength_total = None

for i, year in enumerate(years):

    winter_mask = np.logical_or(
        np.logical_and(ds.month >= 11, ds.year == year),
        np.logical_and(ds.month <= 2, ds.year == year+1)
    )
    
    winter_velocities = ds.v[winter_mask, :, :]
    winter_peak = winter_velocities.max(dim='time')

    middle_mask = np.logical_and(
        np.logical_and(ds.month >= 2.05, ds.month <= 4),
        ds.year == year+1,
    )
    middle_velocities = ds.v[middle_mask, :, :]
    min = middle_velocities.min(dim='time')

    summer_mask = np.logical_and(
        np.logical_and(ds.month >= 4, ds.month <= 6),
        ds.year == year+1,
    )
    summer_velocities = ds.v[summer_mask, :, :]
    summer_peak = summer_velocities.max(dim='time')

    ###### calculate strength of double peaks 
    
    mask = (winter_peak < min) | (summer_peak < 500)
    
    strength = (winter_peak - min)
    #strength = np.where(mask, 0, strength)
    
    if strength_total is None:
        strength_total = strength
    else:
        strength_total += strength

# Average strength over the number of years
strength_avg = strength_total / len(years)

########################
# get strength along centerline

#initialize arrays
strength_centerline = np.zeros((len(x)))

# find indices along centerline 
for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(xx[:, :] - x[i]).argmin(axis = 1)
    target_y_idx= np.abs(yy[:, :] - y[i]).argmin(axis = 0)
    
     # Extract the slope value closest to the target coordinates
    strength_centerline[i] = strength_avg[target_y_idx[0], target_x_idx[0]]
 
########################
# make figures 

fs = 24 #font size

fig, ax3 = plt.subplots(figsize=(16, 16))

# Plot the first dataset
ax3.plot(d, strength_centerline, 'b-')
ax3.set_xlabel('Distance from terminus', fontsize=fs)
ax3.set_ylabel('Strength of winter peak', fontsize=fs, color='b')
ax3.tick_params(axis='y', labelcolor='b', which='major', labelsize=fs)
ax3.tick_params(axis='x', which='major', labelsize=fs)
#second y-axis
ax4 = ax3.twinx()
ax4.plot(d, dem_centerline, 'r-')
ax4.set_ylabel('Elevation [m]', fontsize=fs, color='r')
ax4.tick_params(axis='y', labelcolor='r', which='major', labelsize=fs)
text_objs = plt.gcf().findobj(plt.Text)
for text_obj in text_objs:
    text_obj.set_fontsize(fs)
plt.title('Average from 2018-2021', fontsize=fs)
plt.savefig('avg2018-2021_strength_dem.png')
