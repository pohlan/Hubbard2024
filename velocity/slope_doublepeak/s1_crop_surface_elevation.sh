#!/bin/bash
# 1. Crops the DEM to the geopackage
# 2. Fills a few DEM holes at high elevation.

if test -f "ifsar_hubbardDEM_reproj.tif"; then
    rm ifsar_hubbardDEM_reproj.tif
fi

if test -f "ifsar_hubbardDEM_cropped.tif"; then
    rm ifsar_hubbardDEM_cropped.tif
fi

#gdalwarp -t_srs ESRI:102247 -cutline hubbard_outline.gpkg -crop_to_cutline ifsar_hubbardDEM.tif ifsar_hubbardDEM_cropped.tif # crop 

gdalwarp -t_srs EPSG:3413 ifsar_hubbardDEM.tif ifsar_hubbardDEM_reproj.tif