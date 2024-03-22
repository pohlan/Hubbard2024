#!/bin/bash

if test -f "taku.geo"; then
    rm taku.geo
fi 

eval "$(conda shell.bash hook)"
conda activate glaciome # currently using glaciome environment
python Contour2geo.py -r 500 -i outline/taku_outline.shp -o taku.geo # convert shape file into geo file
conda deactivate

gmsh taku.geo -2

# To go serial, use this:
# ElmerGrid 14 2 mendenhall.msh -autoclean
# ElmerGrid 14 5 mendenhall.msh -autoclean

# To go parallel, use this:
ElmerGrid 14 2 taku.msh -autoclean -metis 4 0 # change first number to determine number of CPUs
ElmerGrid 14 5 taku.msh -autoclean -metis 4 0


