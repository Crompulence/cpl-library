#!/bin/bash
OPEN_FOAM_CASE=openfoam/
LAMMPS_CASE=lammps/

#Clean CFD
cd ${OPEN_FOAM_CASE}
rm -f ../log.openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

#Clean MD
cd ${LAMMPS_CASE}
rm -f ../log.lammps vmd_out.dcd thermo_output.txt particle_dump/*
cd ../

#Run simulation
MD_EXE="$(which lmp_cpl)"
CFD_EXE="$(which CPLSediFOAM)"
cplexec -c 1 "${CFD_EXE} -case ${OPEN_FOAM_CASE} -parallel > log.openfoam" -m 1 "${MD_EXE} < lammps/single.in"
