from osgeo import gdal, ogr
import os

# Input shapefile and raster file paths
shapefile_path = "ifsar_outline.shp"
raster_file_path = "ifsar_hubbardDEM_reproj.tif"
output_raster_path = "ifsar_hubbardDEM_crop.tif"

# Open the shapefile
shapefile = ogr.Open(shapefile_path)
layer = shapefile.GetLayer()

# Get the geometry from the shapefile
feature = layer.GetNextFeature()
geometry = feature.GetGeometryRef()

# Open the raster file
raster = gdal.Open(raster_file_path)

# Get the extent of the geometry
xmin, xmax, ymin, ymax = geometry.GetEnvelope()

# Define the output raster resolution and geotransform
x_res = raster.GetGeoTransform()[1]
y_res = -raster.GetGeoTransform()[5]
geotransform = (xmin, x_res, 0, ymax, 0, -y_res)

# Create a new raster in memory
mem_driver = gdal.GetDriverByName('MEM')
clipped_raster = mem_driver.Create('', int((xmax-xmin)/x_res), 
int((ymax-ymin)/y_res), 1, gdal.GDT_Float32)

# Set the geotransform and projection
clipped_raster.SetGeoTransform(geotransform)
clipped_raster.SetProjection(raster.GetProjection())

# Perform the clip
gdal.RasterizeLayer(clipped_raster, [1], layer, burn_values=[1])

# Create output raster file
driver = gdal.GetDriverByName('GTiff')
out_raster = driver.CreateCopy(output_raster_path, clipped_raster)

# Clean up
del shapefile, raster, clipped_raster, out_raster


