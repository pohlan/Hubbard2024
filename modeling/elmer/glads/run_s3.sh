#!/bin/bash
# Run s1_diagnostic in parallel on four processors.

mpirun -np 4 -quiet ElmerSolver s3_coulomb.sif
