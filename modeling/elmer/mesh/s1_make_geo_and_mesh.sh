#!/bin/bash

if test -f "hubbard.geo"; then
    rm hubbard.geo
fi 

#python Contour2geo.py -r 500 -i outline/hubbard_contour.shp -o hubbard.geo # convert shape file into geo file

python Contour2geo.py -r 150.0 -i ./outline/hubbard_contour.dat -o hubbard_mesh.geo

gmsh hubbard.geo -2 