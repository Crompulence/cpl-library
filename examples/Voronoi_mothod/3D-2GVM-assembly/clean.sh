#!/bin/bash

# Clean CPL
cd cpl
rm -f coupler_header map_CFD map_MD
cd ../

# Clean CFD
cd openfoam/
rm -rf ../log.openfoam processor*/ constant/polyMesh/{boundary,faces,neighbour,owner,points}
rm -rf [0-9].* [1-9]
cd ../

#Clean MD
cd lammps/
rm -f ../log.lammps print_*.txt
rm -rf assembly*
rm -rf *vtk
cd ../

# Clean results
rm -rf __pycache__/ .cache/

# Clean pov files
rm -rf *_p*.pov *_v*.pov 
rm -rf *.zip
cd pov
rm -rf *_p*.pov *_v*.pov
rm -rf *.png
