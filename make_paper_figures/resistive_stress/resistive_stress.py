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
from scipy.integrate import cumulative_trapezoid

# parameters

g = 9.81
rhoi = 917
rhow = 1000

# calculate resulting flow force from velocity cube
A = 24e-25  # (PoG table)
ng = 3  # glens exponent


# load centerline

points = pd.read_csv("../data/centerline_points_100m.csv")
x, y = points.X.to_numpy(), points.Y.to_numpy()

# Reverse the arrays using slicing
x = x[::-1]
y = y[::-1]

# distance from terminus along centerline
d = np.linspace(0, 100 * len(x), len(x))

# spatial resolution of about 1000m (2xice average thickness)
l = np.linspace(0, 100 * len(d), 20)
l_bar = np.repeat(l, 2)

print("finished loading centerlines")


# load DEM
# Open the GeoTIFF file
with rio.open("../data/ifsar_hubbardDEM_reproj.tif") as src:
    # Convert centerline coordinates to pixel indices
    i, j = rio.transform.rowcol(src.transform, x, y)

    # Read the raster data
    data = src.read(1)  # Assuming a single band image

    surfdem_centerline = data[i, j]

print("finished calculating dem and slope along centerline")


# load DEM

with rio.open("../data/hubbard_bedrock_icebridge_reproj.tif") as src:
    # Convert centerline coordinates to pixel indices
    i, j = rio.transform.rowcol(src.transform, x, y)

    # Read the raster data
    data = src.read(1)  # Assuming a single band image

    beddem_centerline = data[i, j]


n = 5
beddem_centerline = np.convolve(np.ones(n) / n, beddem_centerline, mode="same")
beddem_centerline[0:5] = beddem_centerline[6]
beddem_centerline[-5:] = beddem_centerline[-6]

surfdem_centerline = np.convolve(np.ones(n) / n, surfdem_centerline, mode="same")
surfdem_centerline[0:5] = surfdem_centerline[6]
surfdem_centerline[-5:] = surfdem_centerline[-6]

H = surfdem_centerline - beddem_centerline
dhdx = np.gradient(surfdem_centerline, d)
h = surfdem_centerline

# calculate force acting on terminus
Hf = H[0]
hf = h[0]
Ff = g * rhoi / 2 * ((1 - rhow / rhoi) * Hf**2 + rhow / rhoi * hf * (2 * Hf - hf))

print("finished calculating thickness along centerline")


# velocity along centerline at a particular time

hubv = xarray.open_dataset("../Hubbard_5eminus5.nc")
# hubv = xarray.open_dataset("../Hubbard_sentinel1.nc")
year = 2016
ns_in_day = 60 * 60 * 24 * 1e9
epoch = np.datetime64("%s-01-01" % year)
t = ((hubv.time[:] - epoch).to_numpy() / ns_in_day).astype(np.float32)

# Convert time values to datetime objects
datetime_index = pd.to_datetime(hubv.time.values)

# time series for each point
vx = np.zeros((x.shape[0], hubv.vx.shape[0]))
vy = np.zeros((x.shape[0], hubv.vx.shape[0]))
v = np.zeros((x.shape[0], hubv.vx.shape[0]))


# plot different times during the year

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5), sharex=True, tight_layout=True)
fs = 11

# Plot the first set of data on ax1
ax1.plot(d, surfdem_centerline, label="Surface Elevation")
ax1.plot(d, beddem_centerline, label="Bed Elevation", color="black")
ax1.set_ylabel("Elevation [m]", fontsize=fs)
ax1.tick_params(axis="both", which="major", labelsize=fs)


# m = 1000  # September 2018
# m = 1100 # Jan 2019
# m = 1225 # May 2019
# m = 1250 # June 2019
# m =
slices = [1100, 1225, 1350]

lines_2 = []
labels_2 = []
colors = ["yellow", "orange", "red", "purple"]

for j, m in enumerate(slices):
    print(j, m)
    date = datetime_index[m]
    # time series for each point
    vx = np.zeros((x.shape[0], hubv.vx.shape[0]))
    vy = np.zeros((x.shape[0], hubv.vx.shape[0]))
    v = np.zeros((x.shape[0], hubv.vx.shape[0]))

    for i in range(len(x)):
        # get indices of coordinates closest to points of interest
        target_x_idx = np.abs(hubv.x.values - x[i]).argmin()
        target_y_idx = np.abs(hubv.y.values - y[i]).argmin()

        # Extract the time series closest to the target coordinates
        vx[i, :] = hubv.vx[m, target_y_idx, target_x_idx]
        vy[i, :] = hubv.vy[m, target_y_idx, target_x_idx]
        v[i] = np.sqrt(vx[i] ** 2 + vy[i] ** 2)

    vx = vx[:, 0]
    vy = vy[:, 0]
    v = v[:, 0]

    # smooth velocity
    n = 10  # smoothing value
    vy = np.convolve(np.ones(n) / n, vy, mode="same")
    vx = np.convolve(np.ones(n) / n, vx, mode="same")

    vy[0:5] = vy[6]
    vx[0:5] = vx[6]
    vy[-5:] = vy[-6]
    vx[-5:] = vx[-6]
    v_smooth = np.sqrt(vx**2 + vy**2)

    # Compute the gradient for each component
    dvx_dx = np.gradient(vx, d)
    dvy_dx = np.gradient(vy, d)

    # Compute the magnitude of the gradient along centerline
    dvdx = np.sqrt(dvx_dx**2 + dvy_dx**2)
    dvdx = np.convolve(np.ones(n) / n, dvdx, mode="same")

    dvdx[0:5] = dvdx[6]
    dvdx[-5:] = dvdx[-6]

    print("Velocity on", date)

    epsilon0 = dvdx / 2
    eta = (2 * A * ((epsilon0 / A) ** ((ng - 1) / ng))) ** (-1)

    # driving force
    taud = g * rhoi * H * dhdx  # driving stress
    Fd = cumulative_trapezoid(taud, d, initial=0) + Ff

    # resulting longitudinal force
    net_f = 4 * H * dvdx * eta

    # cumulative restisting force
    Rf = Fd - net_f

    # calculate average resisting force
    Rf_avg = np.zeros(len(l))

    for i in range(len(l) - 1):
        Rf_avg[i] = Rf[i + 1] - Rf[i]

    # calculate average driving force
    Fd_avg = np.zeros(len(l))

    for i in range(len(l) - 1):
        Fd_avg[i] = Fd[i + 1] - Fd[i]

    Fd_avg_bar = np.repeat(Fd_avg, 2)
    Rf_avg_bar = np.repeat(Rf_avg, 2)

    # Plot the second set of data on ax2
    (line,) = ax2.plot(
        l_bar[1:],
        Rf_avg_bar[:-1],
        label=f"{date.strftime('%b %Y')}",
        color=colors[j],
    )

    lines_2.append(line)

ax2.set_xlabel("Distance from Terminus", fontsize=fs)
ax2.set_ylabel("Resistive force per unit width", fontsize=fs)
ax2.tick_params(axis="both", which="major", labelsize=fs)
# ax2.set_ylim(np.min(Rf_avg)-1e9, np.max(Rf_avg)+1e9)  # Example limits; adjust as needed

# Adding titles and grid
# plt.suptitle("Elevation and Resisting Force in 2019", fontsize=fs)
ax1.grid(True)
ax2.grid(True)

# Adding legends
lines_1, labels_1 = ax1.get_legend_handles_labels()
ax1.legend(lines_1, labels_1, loc="best")
ax2.legend(loc="best")
# ax2.legend(lines_2, labels_2, loc="best")
plt.savefig("elevation_resist_f.png")
plt.show()
