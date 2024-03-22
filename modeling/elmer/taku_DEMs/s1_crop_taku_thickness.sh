#!/bin/bash
# Converts thickness map to UTM zone 8N, crops to Mendenhall region, and assigns a minimum value of 1 m.

xmin=522500
xmax=561000
ymin=6473500
ymax=6524500

if test -f "taku_thickness.tif"; then
    rm taku_thickness.tif
fi 
 
gdalwarp -t_srs epsg:32608 -te $xmin $ymin $xmax $ymax -tr 50 50 ./Millan/RGI-2/THICKNESS_RGI-2.1_2021July09.tif temp.tif     # crop and reproject
gdal_translate -a_nodata none temp.tif taku_thickness.tif # convert nodata into data

rm temp.tif

