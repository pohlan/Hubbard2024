import numpy as np
import matplotlib.pyplot as plt
import rasterio
import pandas as pd

# open bed DEM 
file_name = '/Users/amyjenson/Documents/Github/Hubbard2024/modeling/elmer/DEMs/hubbard_bed.tif'
with rasterio.open(file_name) as src:
    zb = src.read(1)
    height = zb.shape[0]
    width = zb.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xd, yd = rasterio.transform.xy(src.transform, rows, cols)
    lons_bed= np.array(xd)
    lats_bed = np.array(yd)

xb = lons_bed.flatten()
yb = lats_bed.flatten()

zb = zb.flatten()
zb = np.nan_to_num(zb, nan=-9999)
value = -10000
mat = np.where(zb > value, True, False)
bedDEM = np.where(mat, zb, -9999)
bedDEM = bedDEM.flatten()

# Sort indices based on x and y coordinates
sort_indices = np.lexsort((yb, xb))
sorted_xb = xb[sort_indices]
sorted_yb = yb[sort_indices]
sorted_zb = zb[sort_indices]

#Save the bedace in an ascii file 
np.savetxt('./hubbard_bed.dat',np.column_stack((sorted_xb, sorted_yb, bedDEM)),fmt="%10.4f %10.4f %10.4f")
 

# Calculate parameters for SIF file
x0 = sorted_xb[0]  # bottom left x coordinate
y0 = sorted_yb[0]  # bottom left y coordinate
lx = np.max(sorted_xb) - np.min(sorted_xb)  # the x length of the covered domain
ly = np.max(sorted_yb) - np.min(sorted_yb)  # the y length of the covered domain
nx = len(np.unique(sorted_xb))  # number of levels in x direction
ny = len(np.unique(sorted_yb))  # number of levels in y direction


print('bed, x0 =', f'{x0:.4f}', 'y0 =', f'{y0:.4f}', 'lx =', f'{lx:.4f}', 'ly =', f'{ly:.4f}', 'nx =', nx, 'ny =', ny)

# Read the mesh_parameters.IN file
with open('mesh_parameters.IN', 'r') as file:
    lines = file.readlines()

# Find and replace existing parameters
parameters_to_replace = {
    '# bed_x0': x0,
    '# bed_y0': y0,
    '# bed_lx': lx,
    '# bed_ly': ly,
    '# bed_nx': nx,
    '# bed_ny': ny
}

for i in range(len(lines)):
    for parameter, value in parameters_to_replace.items():
        if parameter in lines[i]:
            lines[i] = f'{parameter} = {value:.4f}\n'
            parameters_to_replace.pop(parameter)  # Remove the parameter once it's replaced
            break  # Move to the next line

# Append new parameters if they were not found in the file
for parameter, value in parameters_to_replace.items():
    lines.append(f'{parameter} = {value:.4f}\n')

# Write the edited content back to mesh_parameters.IN file
with open('mesh_parameters.IN', 'w') as file:
    file.writelines(lines)
