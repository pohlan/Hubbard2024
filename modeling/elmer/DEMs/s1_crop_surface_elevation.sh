#!/bin/bash
# 1. Crops the DEM to the geopackage
# 2. Fills a few DEM holes at high elevation.

if test -f "hubbard_surface.tif"; then
    rm hubbard_surface.tif
fi

gdalwarp -t_srs ESRI:102247 -tr 50 50 -cutline hubbard_outline.gpkg -crop_to_cutline ifsar_hubbardDEM.tif ifsar_hubbardDEM_cropped.tif # crop 

gdal_fillnodata.py ifsar_hubbardDEM_cropped.tif ./hubbard_surface.tif # fill data 

rm ifsar_hubbardDEM_cropped.tif
