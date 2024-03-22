#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from osgeo import gdal
import numpy as np
import pandas as pd
import datetime as dt

from matplotlib import pyplot as plt

from scipy.ndimage import gaussian_filter

#%% correct DEM
file = 'hubbard_thickness.tif'

ds = gdal.Open(file)  # ds = data source
outFileName = './hubbard_thickness_smoothed.tif'

# grab data, determine size of array
data = ds.ReadAsArray()
rows, cols = data.shape

# smooth the bed with gaussian filter
data1 = gaussian_filter(data,5,order=0)



driver = gdal.GetDriverByName("GTiff")
# create new file
outdata = driver.Create(outFileName, cols, rows, 1, gdal.GDT_Float32)
outdata.SetGeoTransform(ds.GetGeoTransform()) # sets same geotransform as input
outdata.SetProjection(ds.GetProjection())  # sets same projection as input
outdata.GetRasterBand(1).WriteArray(data1)
# outdata.GetRasterBand(1).SetNoDataValue()
outdata.FlushCache()  # saves to disk!!

outdata = None
