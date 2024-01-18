#!/bin/bash
gdal_merge.py H130*.bag -o H130_merge.tif
python mask_nodata.py
gdal_translate -a_nodata nan -b 1 -co COMPRESS=LZW H130_merge.tif H130_ellipsoid.tif
rm H130_merge.tif
