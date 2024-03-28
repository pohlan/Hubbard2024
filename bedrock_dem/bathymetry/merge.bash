#!/bin/bash
gdal_merge.py H130*.bag -o H130_merge.tif
python mask_nodata.py
gdal_translate -a_nodata nan -b 1 -co COMPRESS=LZW H130_merge.tif H130_ellipsoid.tif
gdalwarp -tr 50 50 -t_srs EPSG:3338 H130_ellipsoid.tif H130_3338.tif
rm H130_merge.tif
