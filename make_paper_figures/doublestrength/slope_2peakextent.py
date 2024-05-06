import xarray
import matplotlib.pyplot as plt
import numpy as np
import scipy
from tqdm import tqdm
import rasterio as rio
import affine
import os
import elevation
import imageio.v2 as imageio
from matplotlib.colors import BoundaryNorm
import pandas as pd
import geopandas as gpd

# open dataset
year = "2018"
hubv = xarray.open_dataset("hubbard_%s.nc" % year)

# t0 = np.datetime64("%s-09-01" % (str(int(year)-1)))
# t1 = np.datetime64("%s-10-01" % year)
# mask = np.logical_and(hubv.time > t0, hubv.time < t1)
import warnings
ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-01-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)

amp0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
amp1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
phase0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
phase1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))

amp0[:] = np.nan
amp1[:] = np.nan
phase0[:] = np.nan
phase1[:] = np.nan


# Fit smooth function to data and identify peaks
for i in tqdm(range(hubv.speed.shape[1])):
    for j in range(hubv.speed.shape[2]):
        v = hubv.speed[:,i,j]
        ti = np.linspace(np.min(t), np.max(t), 1000)
        
        # Suppress all rank warnings
        warnings.filterwarnings('ignore', category=np.RankWarning)

        # Polynomial fit
        m = np.polyfit(t, v, 12)

        vi = np.polyval(m, ti)

        ##### this looks for mimima instead of maxima
        #vi = -vi 
        pk = scipy.signal.find_peaks(vi, prominence=50)[0]
        #vi = -vi 
        ######
        
        if(len(pk) == 1):
            amp0[i,j] = vi[pk[0]]
            phase0[i,j] = ti[pk[0]]
            continue
        elif(len(pk) == 2):
            amp0[i,j] = vi[pk[0]]
            phase0[i,j] = ti[pk[0]]
            amp1[i,j] = vi[pk[1]]
            phase1[i,j] = ti[pk[1]]
            continue


ds = xarray.load_dataset("hubbard_%s.nc" % year)

res = ds.x[1]-ds.x[0]

xform = affine.Affine.translation(ds.x[0] - res / 2, ds.y[0] - res / 2) * affine.Affine.scale(res, -res)

#v = np.sqrt(ds.vx**2 + ds.vy**2)
#vmean = np.mean(v, axis=0)

with rio.open(
    "./amp1.tif",
    "w",
    driver='GTiff',
    height=amp1.shape[0],
    width=amp1.shape[1],
    count=1,
    dtype=amp1.dtype,
    crs="EPSG:3413",
    transform=xform,
) as ds_out:
    ds_out.write(amp1, 1)


# get extent of double peak plot

# Given center coordinate
center_x = -3310000
center_y = 257700

# Offset for moving northward (in meters)
offset = 8000

xmin = center_x - offset
xmax = center_x + offset
ymin = center_y - offset
ymax = center_y + offset



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

#####
# read file and calculate slope
image = imageio.imread('ifsar_hubbardDEM_reproj.tif')
slope = calculate_slope_angle(image, 5)

df = pd.read_csv("peak_outline.csv") 
df["map_x"] = df["X"]  # Replace "x_column_name" with the actual name of your x-column
df["map_y"] = df["Y"]  # Replace "y_column_name" with the actual name of your y-column


# Adjust the figure size by specifying figsize=(width, height)
fig, axs = plt.subplots(1, 1, figsize=(12, 6))  # Increase width and height as needed

im0 = axs.imshow(slope, extent=[left, right, bottom, top], vmin=0, vmax=20, cmap="viridis")
axs.set_xlim([xmin, xmax])
axs.set_ylim([ymin, ymax])
axs.set_xticks([])
axs.set_yticks([])
axs.plot(df["map_x"], df["map_y"], color='r')

# Set the title of the entire figure
fig.suptitle("Slope and extent of double peak in %s " % year)

# Add colorbars to the subplots
fig.colorbar(im0, ax=axs, label="Slope (deg)")

# Save the figure
plt.savefig("%s_Hubbard_Slope.pdf" % year, bbox_inches="tight")

