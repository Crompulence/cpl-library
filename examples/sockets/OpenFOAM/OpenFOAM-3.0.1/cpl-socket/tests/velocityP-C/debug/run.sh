#mpirun -np 27 python dummyMD.py : -np 27 CPLIcoFoam -parallel -case openfoam
mpirun -n 27 python dummyMD_vels.py & PID=$!; mpirun -n 27 CPLIcoFoam -parallel -case test_vels_case ; wait $PID
reconstructPar -case test_vels_case
stressComponents -case test_vels_case
writeCellCentres -case test_vels_case
