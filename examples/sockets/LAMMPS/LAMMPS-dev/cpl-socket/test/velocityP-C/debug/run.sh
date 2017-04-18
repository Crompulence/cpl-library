#mpiexec -n 27 lmp_cpl -in lammps_forces.in : -n 27 python dummyCFD_forces.py
mpiexec -n 27 lmp_cpl -in lammps_vels.in & PID=$!; mpiexec -n 27 python dummyCFD_vels.py ; wait $PID
