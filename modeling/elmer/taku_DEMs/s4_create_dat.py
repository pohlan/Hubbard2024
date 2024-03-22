# -*- coding: utf-8 -*-
# Create a geo (gmsh input file) file from a contour file
# the contour file contains the (x,y) coordinates of the ordered
# points defining the contour of the domain

import numpy as np
from osgeo import gdal

#%% converts geotiffs into dat files for loading into Elmer/Ice
def createDat(file):
    ds = gdal.Open(file)
    gt = ds.GetGeoTransform()

    elevation = np.flipud(ds.ReadAsArray())

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
    elev_flat = elevation.flatten()

    np.savetxt(file[:-4] + '.dat',np.column_stack((x_flat, y_flat, elev_flat)), fmt="%10.2f %10.2f %10.2f")

    return(west, east, north, south, )

#%%
# create dat files

bedDEM = './taku_bed.tif'
createDat(bedDEM)

surfDEM = './taku_surface.tif'
west, east, north, south = createDat(surfDEM)

