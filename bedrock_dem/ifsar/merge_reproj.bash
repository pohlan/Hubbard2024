#!/bin/bash
gdal_merge.py -o tmp0.tif USGS_AK5M_Alaska_Mid_Accuracy_DEM_Summer_2015_DTM_N6000W13915P.tif USGS_AK5M_Alaska_Mid_Accuracy_DEM_Summer_2015_DTM_N6000W13930P.tif USGS_AK5M_Alaska_Mid_Accuracy_DEM_Summer_2015_DTM_N6000W13945P.tif USGS_AK5M_Upper_Southeast_Alaska_Mid_Accuracy_DEM_1684.tif
gdalwarp -tr 50 50 -t_srs EPSG:3338 tmp0.tif USGS_AK5M_Alaska_Mid_Accuracy_DEM_Summer_2015_DTM_Hubbard.tif
rm tmp0.tif
