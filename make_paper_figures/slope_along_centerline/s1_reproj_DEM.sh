#!/bin/bash
# 1. Reproject 

if test -f "ifsar_hubbardDEM_reproj.tif"; then
    rm ifsar_hubbardDEM_reproj.tif
fi

gdalwarp -t_srs EPSG:3413 -tr 5 5 ifsar_hubbardDEM.tif ifsar_hubbardDEM_reproj.tif