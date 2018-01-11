mpiexec -n 9 python2 md_test.py : -n 9 python2 cfd_test.py
#mpiexec -n 27 lmp_cpl -in lammps_vels.in & PID=$!; mpiexec -n 27 python2 dummyCFD_vels.py ; wait $PID
