#!/bin/bash

# Clean directory
./clean.sh

# Build CFD grid
cd openfoam/
blockMesh
decomposePar
cd ../

#Run simulation
cplexec -c 1 "CPLSediFOAM -case openfoam/ -parallel > log.openfoam" -m 1 "lmp_cpl < lammps/column.in" 
# cplexec -c 1 "CPLCFDDEMFoam -case openfoam/ -parallel > log.openfoam" -m 1 "lmp_cpl < lammps/column.in" 