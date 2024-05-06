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
from numpy.polynomial.polynomial import Polynomial

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))


########################
# open dataset
year = "2018"
hubv = xarray.open_dataset("hubbard_%s.nc" % year)
x_coords = hubv.x
y_coords = hubv.y

########################
# load centerline

points = pd.read_csv('centerline_points_100m.csv')

#transformer = Transformer.from_crs("epsg:3413", "epsg:32607") # UTM 6N
#x, y = transformer.transform(points.X.to_numpy(),points.Y.to_numpy())

x, y =points.X.to_numpy(), points.Y.to_numpy()

terminus = (np.min(x), np.max(y))

#distance from terminus along centerline 
d = np.linspace(100*len(x), 0, len(x))

target_x_idx = np.zeros(len(x), dtype=int)
target_y_idx = np.zeros(len(x), dtype=int)

print('finished loading centerlines')
    
########################
# get strength of double peak along centerline

ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-01-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)

# Convert time values to datetime objects
datetime_index = pd.to_datetime(hubv.time.values)

# time series for each point 
time_series = np.zeros((len(x), hubv.speed.shape[0]))

double_strength = np.zeros((len(x)), dtype=float)

for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    target_x_idx[i] = np.abs(x_coords - x[i]).argmin()
    target_y_idx[i]= np.abs(y_coords - y[i]).argmin()
    
    # Extract the time series closest to the target coordinates
    time_series[i, :] = hubv.speed[:, target_y_idx[i], target_x_idx[i]]
    
    # fit a polynomial of order 4
# this shows the full polynomial, 
    m_poly = Polynomial.fit(t, time_series[i,:], 4) # shape of m is (5,)
    
    # this is an array of coefficients with the form [a0, a1, a2, a3, a4] where a0 + a1 x + a2 x^2 + a3 x^3  + a4 x^4
    m = m_poly.coef
    #print(m)
    
    # Check if m[2] > 0 and m[4] < 0
    if m[2] > 0 and m[4] < 0:
        # Calculate double_strength using the provided formula
        double_strength[i] = (m[2] / np.abs( m[4]))  - (np.abs(m[1])+np.abs(m[3]))
    else:
        # If the condition is not satisfied, set double_strength to zero
        double_strength[i] = 0
 
 #### other ways
 
    # store values of a2 which should give strength of double peak 
    #double_strength[i] = m[2]
    #a2 and 1/a4 which should give strength of double peak 
    double_strength[i] = (m[2]) + (-1/m[4])
 
#double_strength[double_strength == -np.inf] = np.nan    
#double_strength[double_strength == np.inf] = np.nan    

# Split the array into positive and negative parts
#positive_values = double_strength[double_strength > 0]
#negative_values = double_strength[double_strength < 0]

# Scale positive values to fit between 0 and 1
#scaled_positive_values = positive_values / np.max(positive_values)

# Scale negative values to fit between -1 and 0
#scaled_negative_values = negative_values / np.abs(np.min(negative_values))

# Combine scaled positive and negative values
#scaled_double_strength = np.zeros_like(double_strength, dtype=float)
#scaled_double_strength[double_strength > 0] = scaled_positive_values
#scaled_double_strength[double_strength < 0] = scaled_negative_values

#scaled_double_strength = (double_strength - np.nanmin(double_strength))/ (np.nanmax(double_strength) - np.nanmin(double_strength))

scaled_double_strength = (double_strength)/ (np.max(double_strength))

######

print('finished calculating strength of double peak along centerline')

#zero = np.zeros_like(d) 

########################

ax.plot(d, scaled_double_strength)

#ax.plot(d, zero)

# Set x-axis label
ax.set_xlabel('Distance from terminus', fontsize=fs)

# Set y-axis label
ax.set_ylabel('Strength of double peak', fontsize=fs)

#ax.invert_xaxis()

# Set font size for tick labels
ax.tick_params(axis='both', which='major', labelsize=fs)

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('doublepeak_strength.png')