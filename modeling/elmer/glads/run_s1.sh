#!/bin/bash
# Run s1_diagnostic in parallel on four processors.
#mpirun -np 4 -quiet ElmerSolver s1_initial_diagnostic.sif

# Run in serial 
ElmerSolver s1_initial_diagnostic.sif