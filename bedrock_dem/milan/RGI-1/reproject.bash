#!/bin/bash
gdalwarp -r bilinear -tr 50 50 -t_srs EPSG:3338 THICKNESS_RGI-1.1_2021July01.tif THICKNESS_RGI-1.1_2021July01_3338.tif
