#!/usr/bin/env python
from mpi4py import MPI
from cplpy import CPL

comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([2, 2, 2]), xyzL=[1.0, 1.0, 1.0], 
             xyz_orig=[0.0, 0.0, 0.0])
recv_array, send_array = CPL.get_arrays(recv_size=1, send_size=4)

for time in range(5):

    send_array[0,:,:,:] = 5.*time
    CPL.send(send_array)
    recv_array, ierr = CPL.recv(recv_array)
    print("MD", time, recv_array[0,0,0,0])

CPL.finalize()
MPI.Finalize()
