#!/bin/bash

if test -f "hubbard.msh"; then
    rm hubbard.msh
fi 

gmsh hubbard.geo -2 