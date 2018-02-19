import numpy as np
from mpi4py import MPI
from cplpy import CPL

g = 9.81

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)

## Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([1.5E-003, 1.5E-003, 2.5E-003], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([8, 8, 8], order='F', dtype=np.int32)

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=9)

ft = True
for time in range(101):

    # send data to update
    send_array[2,:,:,:] = -5.9490638385009208e-08*g  # * mi
    CPL.send(send_array)

    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array)

    print(time)

CPL.finalize()
MPI.Finalize()




