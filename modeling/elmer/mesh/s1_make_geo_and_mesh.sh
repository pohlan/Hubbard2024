#!/bin/bash

if test -f "hubbard.geo"; then
    rm hubbard.geo
fi 

#python Contour2geo.py -r 500 -i outline/hubbard_ice_contour.shp -o hubbard_ice.geo # convert shape file into geo file

python Contour2geo.py -r 150.0 -i ./outline/hubbard_ice_contour.dat -o hubbard_ice_mesh.geo

gmsh hubbard_ice_mesh.geo -2 