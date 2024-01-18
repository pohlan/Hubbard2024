from osgeo import gdal
import numpy as np

path = "./H130_merge.tif"
ds = gdal.Open(path, gdal.GA_Update)

b1 = ds.GetRasterBand(1).ReadAsArray()
b2 = ds.GetRasterBand(2).ReadAsArray()
mask = np.logical_and(b1==0, b2==0)

b1[mask] = np.nan

ds.GetRasterBand(1).WriteArray(b1)

del ds