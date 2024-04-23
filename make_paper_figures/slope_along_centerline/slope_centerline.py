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

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))


########################
# open dataset
year = "2018"
hubv = xarray.open_dataset("hubbard_%s.nc" % year)


########################
# load centerline

points = pd.read_csv('centerline_points_100m.csv')
transformer = Transformer.from_crs("epsg:3413", "epsg:32607") # UTM 6N
x, y = transformer.transform(points.X.to_numpy(),points.Y.to_numpy())


########################
# get slope along centerline 


def calculate_slope_angle(dem, cell_size):
    """
    Calculate slope angle from a Digital Elevation Model (DEM).
    
    Parameters:
        dem (numpy.ndarray): 2D array representing the DEM.
        cell_size (float): Size of each cell in the DEM (in meters).
        
    Returns:
        numpy.ndarray: 2D array containing the slope angles (in degrees).
    """
    # Calculate gradient using central differences
    dz_dx, dz_dy = np.gradient(dem, cell_size, cell_size)
    
    # Calculate slope angle
    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    slope_deg = np.degrees(slope_rad)
    
    return slope_deg
    
    
# Open the GeoTIFF file
with rio.open('ifsar_hubbardDEM_reproj.tif') as src:
    # Read the raster data
    data = src.read(1)  # Assuming a single band image
    
 # Get the affine transformation matrix
    transform = src.transform
    left, bottom, right, top = src.bounds

# Generate x and y coordinates for each pixel
    rows, cols = data.shape
    x_coords, y_coords = [], []
    for r in range(rows):
        for c in range(cols):
            x, y = rio.transform.xy(transform, r, c)
            x_coords.append(x)
            y_coords.append(y)

    # Store DEM values
    dem = data.flatten()


#####
# read file and calculate slope
#image = imageio.imread('ifsar_hubbardDEM_reproj.tif')
slope = calculate_slope_angle(dem, 5)

# find indices along centerline 
slope_centerline = np.zeros((x.shape[0], ))

for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(x_coords - x[i]).argmin()
    target_y_idx = np.abs(y_coords - y[i]).argmin()
    
    # Extract the time series closest to the target coordinates
    slope_centerline[i, :] = slope[:, target_y_idx, target_x_idx]



########################
# get strength of double peak along centerline

ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-01-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)

# Convert time values to datetime objects
datetime_index = pd.to_datetime(hubv.time.values)

# time series for each point 
time_series = np.zeros((x.shape[0], hubv.speed.shape[0]))

for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(hubv.x.values - x[i]).argmin()
    target_y_idx = np.abs(hubv.y.values - y[i]).argmin()
    
    # Extract the time series closest to the target coordinates
    time_series[i, :] = hubv.speed[:, target_y_idx, target_x_idx]





########################

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

plt.savefig('centerline_velocities.png')