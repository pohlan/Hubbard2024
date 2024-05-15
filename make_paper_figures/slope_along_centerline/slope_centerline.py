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

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))

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

ds = xarray.load_dataset("../Hubbard_5eminus5.nc")
ds["v"] = np.sqrt(ds.vx**2 + ds.vy**2).fillna(0)
ds["month"] = ds.time.dt.month
ds["year"] = ds.time.dt.year
ds["doy"] = ds.time.dt.dayofyear

for i in range(len(points_X)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(hubv.x.values - points_X[i]).argmin()
    target_y_idx = np.abs(hubv.y.values - points_Y[i]).argmin()
    
    # Extract the time series closest to the target coordinates
    vx[i, :] = hubv.vx[:, target_y_idx, target_x_idx]
    vy[i, :] = hubv.vy[:, target_y_idx, target_x_idx]
    v[i] = np.sqrt(vx[i]**2 + vy[i]**2)


for i, year in enumerate([2018]):

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
    strength = np.where(mask, 0, strength)
 
########################
# get slope along centerline

#initialize arrays
slope_centerline = np.zeros((len(x)))
dem_centerline = np.zeros((len(x)))
strength_centerline = np.zeros((len(x)))


# find indices along centerline 
for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(x_coords[:, :] - x[i]).argmin(axis = 1)
    target_y_idx= np.abs(y_coords[:, :] - y[i]).argmin(axis = 0)
    
     # Extract the slope value closest to the target coordinates
    slope_centerline[i] = slope[target_y_idx[0], target_x_idx[0]]
    dem_centerline[i] = data[target_y_idx[0], target_x_idx[0]]
    strength_centerline[i] = strength[target_y_idx[0], target_x_idx[0]]

# Combine target_x_idx, target_y_idx, and slope into a DataFrame
df = pd.DataFrame({'x': x, 'y': y, 'slope': slope_centerline})

# Write the DataFrame to a CSV file
df.to_csv('slope_along_centerline.csv')

print('finished calculating slope along centerline')


########################

 
########################

ax.plot(d, slope_centerline)

# Set x-axis label
ax.set_xlabel('Distance from terminus', fontsize=fs)

# Set y-axis label
ax.set_ylabel('Slope [deg]', fontsize=fs)

# Set font size for tick labels
ax.tick_params(axis='both', which='major', labelsize=fs)

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('slope_centerline.png')



ax.plot(d, dem_centerline)

# Set x-axis label
ax.set_xlabel('Distance from terminus', fontsize=fs)

# Set y-axis label
ax.set_ylabel('Elevation [m]', fontsize=fs)

# Set font size for tick labels
ax.tick_params(axis='both', which='major', labelsize=fs)

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('dem_centerline.png')