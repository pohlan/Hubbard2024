import numpy as np
import matplotlib.pyplot as plt
import rasterio
import pandas as pd


# open surface DEM 
file_name = '/Users/amyjenson/Documents/Github/Hubbard2024/modeling/elmer/DEMs/hubbard_surface.tif'
with rasterio.open(file_name) as src:
    zs = src.read(1)
    height = zs.shape[0]
    width = zs.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xd, yd = rasterio.transform.xy(src.transform, rows, cols)
    lons_surf= np.array(xd)
    lats_surf = np.array(yd)

#xs = np.flip(lons_surf.flatten())
#ys = np.flip(lats_surf.flatten())
#zs = np.flip(zs.flatten())
xs = lons_surf.flatten()
ys = lats_surf.flatten()
zs = zs.flatten()
value = 0
mat = np.where(zs > value, True, False)
surfDEM = np.where(mat, zs, -9999)

# Sort indices based on x and y coordinates
sort_indices = np.lexsort((ys, xs))
sorted_xs = xs[sort_indices]
sorted_ys = ys[sort_indices]
sorted_zs = zs[sort_indices]

#Save the surface in an ascii file 
np.savetxt('./hubbard_surf.dat',np.column_stack((sorted_xs, sorted_ys, surfDEM)),fmt="%10.4f %10.4f %10.4f")
 

# Calculate parameters for SIF file
x0 = sorted_xs[0]  # bottom left x coordinate
y0 = sorted_ys[0]  # bottom left y coordinate
lx = np.max(sorted_xs) - np.min(sorted_xs)  # the x length of the covered domain
ly = np.max(sorted_ys) - np.min(sorted_ys)  # the y length of the covered domain
nx = len(np.unique(sorted_xs))  # number of levels in x direction
ny = len(np.unique(sorted_ys))  # number of levels in y direction


print('surface, x0 =', f'{x0:.4f}', 'y0 =', f'{y0:.4f}', 'lx =', f'{lx:.4f}', 'ly =', f'{ly:.4f}', 'nx =', nx, 'ny =', ny)

# Read the mesh_parameters.IN file
with open('mesh_parameters.IN', 'r') as file:
    lines = file.readlines()

# Find and replace existing parameters
parameters_to_replace = {
    '# surf_x0': x0,
    '# surf_y0': y0,
    '# surf_lx': lx,
    '# surf_ly': ly,
    '# surf_nx': nx,
    '# surf_ny': ny
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
