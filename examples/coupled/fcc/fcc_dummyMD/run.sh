#Set OpenFOAM folder
OPEN_FOAM_CASE=openfoam/
LAMMPS_CASE=lammps/

#Clean CFD
cd ${OPEN_FOAM_CASE}
python clean.py -f
rm -f ../log.openfoam
blockMesh
decomposePar
cd ../

#Run job
cplexec -c 1 "CPLSediFOAM -case ${OPEN_FOAM_CASE} -parallel > log.openfoam" -m 1 "./MD_dummy_fcc.py"

#Plot results
python test_fcc.py