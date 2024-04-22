import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray
import scipy
from tqdm import tqdm

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

# load points of interest 

points_x = np.array([points.X[1], points.X[3], points.X[5] ])   # specify your target x coordinate
points_y = np.array([points.Y[1], points.Y[3], points.Y[5]] ) # specify your target y coordinate

########################
# plug in these coordinates to victors velocity datacube

ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-01-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)
print(t)

# want a time series for each point 
time_series = np.zeros((points_x.shape[0], v.shape[0]))
#print(time_series.shape)

# Fit smooth function to data and identify peaks
for i in tqdm(range(hubv.speed.shape[1])):
    for j in range(hubv.speed.shape[2]):
        v = hubv.speed[:,i,j]
        

# Calculate the Euclidean distance between each (x, y) coordinate and the target coordinates
for i in range(len(points_x)):
    # Convert the target coordinates to the same coordinate system as your velocity data cube
    target_x_idx = np.abs(v.x.values - points_x[i]).argmin()
    target_y_idx = np.abs(v.y.values - points_y[i]).argmin()

    # Find the index of the minimum distance
    distances = np.sqrt((v.x.values - points_x[i])**2 + (v.y.values - points_y[i])**2)

    # Find the index of the minimum distance
    min_index = np.argmin(distances)

    # Extract the time series closest to the target coordinates
    time_series[i, :] = v[:, target_y_idx, target_x_idx]

    #print(time_series) 
    
    # plot time series at points 
    ax.plot(t, time_series[i,:]) 
    plt.show

########################



########################

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('centerline_velocities.png')