#Set OpenFOAM folder
LAMMPS_CASE=lammps/

#Clean MD
cd ${LAMMPS_CASE}
rm -f ../log.lammps vmd_out.dcd fcc*.dump
cd ../

#Run job
cplexec -c 1 "./CFD_dummy_fcc.py" -m 1 "lmp_cpl -in lammps/fcc.in"

#Plot results
python test_fcc.py
