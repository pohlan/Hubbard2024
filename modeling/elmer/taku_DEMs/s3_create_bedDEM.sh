#!/bin/bash
# Subtracts the ice thickness map from the surface DEM to produce a bedrock DEM

if test -f "taku_bed.tif"; then
    rm taku_bed.tif
fi

A=./taku_surface.tif
B=./taku_thickness_smoothed.tif

gdal_calc.py -A $A -B $B --outfile=./taku_bed.tif --calc="A-B"
