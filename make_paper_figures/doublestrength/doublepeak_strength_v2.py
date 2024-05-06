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
import warnings

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))


########################
# open dataset
#hubv = xarray.open_dataset("../Hubbard_5eminus5")

# open dataset
year = "2018"
hubv = xarray.open_dataset("hubbard_%s.nc" % year)

x_coords = hubv.x
y_coords = hubv.y

########################
# load centerline

points = pd.read_csv('centerline_points_100m.csv')

x, y =points.X.to_numpy(), points.Y.to_numpy()

terminus = (np.min(x), np.max(y))

#distance from terminus along centerline 
d = np.linspace(100*len(x), 0, len(x))

idx = np.zeros(len(x), dtype=int)
idy = np.zeros(len(x), dtype=int)

for i in range(len(x)):
    # get indices of coordinates closest to points of interest
    idx[i] = np.abs(x_coords - x[i]).argmin()
    idy[i]= np.abs(y_coords - y[i]).argmin()
    
print('finished loading centerlines')
    
########################
# get strength of double peak along centerline
ns_in_day = 60*60*24*1e9
epoch = np.datetime64("%s-10-01" % year)
t = ((hubv.time[:]-epoch).to_numpy()/ns_in_day).astype(np.float32)

t0 = np.datetime64("%s-10-01" % (str(int(year)-1)))
t1 = np.datetime64("%s-10-01" % year)
mask = np.logical_and(hubv.time > t0, hubv.time < t1)

max_amp0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
max_amp1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
max_phase0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
max_phase1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))

max_amp0[:] = np.nan
max_amp1[:] = np.nan
max_phase0[:] = np.nan
max_phase1[:] = np.nan

min_amp0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
min_phase0 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
min_amp1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
min_phase1 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
min_amp2 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
min_phase2 = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))

min_amp0[:] = np.nan
min_phase0[:] = np.nan
min_amp1[:] = np.nan
min_phase1[:] = np.nan
min_amp2[:] = np.nan
min_phase2[:] = np.nan

# Convert time values to datetime objects
datetime_index = pd.to_datetime(hubv.time.values)

avg_mag = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))
double_strength = np.zeros((hubv.speed.shape[1], hubv.speed.shape[2]))

for i in tqdm(range(hubv.speed.shape[1])):
    for j in range(hubv.speed.shape[2]):    
        v = hubv.speed[:, i, j]
    
        ti = np.linspace(np.min(t), np.max(t), 1000)
        
        # Suppress all rank warnings
        warnings.filterwarnings('ignore', category=np.RankWarning)

        # Polynomial fit
        m = np.polyfit(t, v, 12)
        vi = np.polyval(m, ti)

        # this looks for maxima
        max_pk = scipy.signal.find_peaks(vi, prominence=50)[0]
        
        # this looks for mimima 
        vi = -vi 
        min_pk = scipy.signal.find_peaks(vi, prominence=1, width=200)[0]
        vi = -vi

        if(len(max_pk) == 1):
            max_amp0[i, j] = vi[max_pk[0]]
            max_phase0[i, j] = ti[max_pk[0]]
            continue
        elif(len(max_pk) == 2):
            max_amp0[i, j] = vi[max_pk[0]]
            max_phase0[i, j] = ti[max_pk[0]]
            max_amp1[i, j] = vi[max_pk[1]]
            max_phase1[i, j] = ti[max_pk[1]]
            avg_mag[i,j] = (max_amp0[i,j] + max_amp1[i,j])/2  # average of two double peaks
            
            continue

        if(len(min_pk) == 1):
            min_amp0[i, j] = vi[min_pk[0]]
            min_phase0[i, j] = ti[min_pk[0]]
            continue
        elif(len(min_pk) == 2):
            min_amp0[i, j] = vi[min_pk[0]]
            min_phase0[i, j] = ti[min_pk[0]]
            min_amp1[i, j] = vi[min_pk[1]]
            min_phase1[i, j] = ti[min_pk[1]]
        elif(len(min_pk)==3):
            min_amp0[i, j] = vi[min_pk[0]]
            min_phase0[i, j] = ti[min_pk[0]]
            min_amp1[i, j] = vi[min_pk[1]]
            min_phase1[i, j] = ti[min_pk[1]]
            min_amp2[i, j] = vi[min_pk[2]]
            min_phase2[i, j] = ti[min_pk[2]] 
            double_strength[i,j] =  (avg_mag[i,j] - min_amp1[i,j])/avg_mag[i,j]
            continue

        ######

plt.figure(figsize=(8, 6))
plt.imshow(min_amp1, cmap='terrain', aspect='auto')
plt.colorbar(label='Strength of double peak')
plt.grid(True)
plt.show()
   
print('finished calculating strength of double peak along centerline')


########################

#distance from terminus along centerline 
d = np.linspace(100*len(x), 0, len(double_strength[0]))

#ax.plot(d, double_strength)

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