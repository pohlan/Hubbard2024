import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray
import scipy
#from tqdm import tqdm
from scipy.signal import savgol_filter

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))


########################
# open dataset
year = "2018"
hubv = xarray.open_dataset("hubbard_%s.nc" % year)


########################
# load coordinates of points for velocity plots 
points = pd.read_csv('centerline_points_3000m.csv') # this is in ESPG: 3413 

#coordinates of points of interest 
points_X = np.array([points.X[2], points.X[3], points.X[4], points.X[5], points.X[6] ])   # specify your target x coordinate
points_Y = np.array([points.Y[2], points.Y[3], points.Y[4], points.Y[5], points.Y[6] ]) # specify your target y coordinate

# dense points around ice fall 

# load coordinates of points for velocity plots 
#points = pd.read_csv('centerline_points_1000m.csv') # this is in ESPG: 3413 

#points_X = np.array([points.X[6], points.X[7], points.X[8], points.X[9], points.X[10] ])   # specify your target x coordinate
#points_Y = np.array([points.Y[6], points.Y[7], points.Y[8], points.Y[9], points.Y[10] ]) # specify your target y coordinate


########################
# get velocities at these locations from velocity datacube

ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-10-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)

# Convert time values to datetime objects
datetime_index = pd.to_datetime(hubv.time.values)


# time series for each point 
time_series = np.zeros((points_X.shape[0], hubv.speed.shape[0]))

for i in range(len(points_X)):
    # get indices of coordinates closest to points of interest
    target_x_idx = np.abs(hubv.x.values - points_X[i]).argmin()
    target_y_idx = np.abs(hubv.y.values - points_Y[i]).argmin()
    
    # Extract the time series closest to the target coordinates
    time_series[i, :] = hubv.speed[:, target_y_idx, target_x_idx]

    # plot time series at points
    #ax.plot(datetime_index, time_series[i, :],  linewidth=2)
    
   # Smooth the time series using Savitzky-Golay filter
    smoothed_time_series = savgol_filter(time_series[i, :], window_length=7, polyorder=3)

    # plot smoothed time series at points
    ax.plot(datetime_index, smoothed_time_series, linewidth=4)


# Set x-axis label
ax.set_xlabel('Time', fontsize=fs)

# Set y-axis label
ax.set_ylabel('Velocity [m/a]', fontsize=fs)

# Set font size for tick labels
ax.tick_params(axis='both', which='major', labelsize=fs)

# Format x-axis to display month/year
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%m/%Y'))

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)


########################

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('centerline_velocities.png')