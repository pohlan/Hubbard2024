#!/bin/bash
# 1. Crops the DEM to the geopackage
# 2. Fills a few DEM holes at high elevation.

if test -f "hubbard_bed_icebridge.tif"; then
    rm hubbard_bed_icebridge.tif
fi

#gdalwarp -t_srs ESRI:102247 -tr 50 50 -cutline hubbard_outline.gpkg -crop_to_cutline hubbard_bedrock_icebridge.tif hubbard_bedrock_icebridge_cropped.tif 

#gdal_fillnodata.py hubbard_bedrock_icebridge_cropped.tif ./hubbard_bed_icebridge.tif # fill data 

#rm hubbard_bedrock_icebridge_cropped.tif

gdalwarp -t_srs ESRI:102247 -tr 50 50 -cutline hubbard_outline.gpkg -crop_to_cutline hubbard_bedrock_icebridge.tif hubbard_bed_icebridge.tif 