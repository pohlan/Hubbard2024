# -*- coding: utf-8 -*-

import numpy as np
from osgeo import gdal

#%% converts geotiffs into dat files for loading into Elmer/Ice
def createDat(file):
    ds = gdal.Open(file)
    gt = ds.GetGeoTransform()

    velocity = np.flipud(ds.ReadAsArray())

    west = gt[0]
    dx = gt[1]
    east = west + dx*ds.RasterXSize
    xL = east - west # length in east-west direction

    north = gt[3]
    dy = gt[5]
    south = north + dy*ds.RasterYSize
    yL = north-south # length in north-south direction

    x = np.linspace(west+dx/2, east-dx/2, ds.RasterXSize, endpoint=True)
    y = np.linspace(south-dy/2, north+dy/2, ds.RasterYSize, endpoint=True)
    X, Y = np.meshgrid(x,y)

    x_flat = X.flatten()
    y_flat = Y.flatten()
    vel_flat = velocity.flatten()

    np.savetxt(file[:-4] + '.dat',np.column_stack((x_flat, y_flat, vel_flat)), fmt="%10.2f %10.2f %10.2f")
    
    # Calculate parameters for SIF file
    x0 = x_flat[0]  # bottom left x coordinate
    y0 = y_flat[0]  # bottom left y coordinate
    lx = np.max(x_flat) - np.min(x_flat)  # the x length of the covered domain
    ly = np.max(y_flat) - np.min(y_flat)  # the y length of the covered domain    
    nx = len(np.unique(x_flat))  # number of levels in x direction
    ny = len(np.unique(y_flat))  # number of levels in y direction


    print('surface, x0 =', f'{x0:.4f}', 'y0 =', f'{y0:.4f}', 'lx =', f'{lx:.4f}', 'ly =', f'{ly:.4f}', 'nx =', nx, 'ny =', ny)

# Read the mesh_parameters.IN file
    with open('mesh_parameters.IN', 'r') as file:
        lines = file.readlines()

# Find and replace existing parameters
    parameters_to_replace = {
        '# x0': x0,
        '# y0': y0,
        '# lx': lx,
        '# ly': ly,
        '# nx': nx,
        '# ny': ny
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

    
    return(west, east, north, south, )

#%%
# create dat files

millan_v = './millan_v.tif'
createDat(millan_v)

millan_vx = './millan_vx.tif'
createDat(millan_vx)

millan_vy = './millan_vy.tif'
createDat(millan_vy)
