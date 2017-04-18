#mpiexec -n 27 lmp_cpl -in lammps_forces.in : -n 27 python dummyCFD_forces.py
mpiexec -n 27 lmp_cpl -in lammps_forces.in & PID=$!; mpiexec -n 27 python dummyCFD_forces.py ; wait $PID
