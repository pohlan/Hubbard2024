#!/bin/bash

gdal_calc.py --extent=intersect -A ifsar/USGS_AK5M_Alaska_Mid_Accuracy_DEM_Summer_2015_DTM_Hubbard.tif -B milan/RGI-1/THICKNESS_RGI-1.1_2021July01_3338.tif --outfile=tmp0.tif --calc="A-B"
gdal_merge.py -a_nodata nan -o tmp1.tif tmp0.tif bathymetry/H130_3338.tif
gdalwarp -cutline term_mask.gpkg -crop_to_cutline tmp1.tif tmp2.tif
gdal_fillnodata.py -md 20 -si 5 tmp2.tif hubbard_bedrock.tif
rm tmp0.tif tmp1.tif tmp2.tif
