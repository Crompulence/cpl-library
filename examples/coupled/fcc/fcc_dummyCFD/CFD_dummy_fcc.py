import numpy as np
from mpi4py import MPI
from cplpy import CPL

#initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([2., 10., 2.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([5, 200, 5], order='F', dtype=np.int32)

#initialise CPL
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)
print("After CPL setup")

#Setup send and recv buffers
recv_array, send_array = CPL.get_arrays(recv_size=8, send_size=9)

#Main time loop
porousStart = 91
porousEnd = 108
phi = 0.74
Ub = 0.01;
K = 435.
mu = 1.e-2
dp = 0.0565685

eps = 1. - phi
rp = 0.5*dp
gradP = 4.5*mu*phi*K*Ub/(rp**2)

for time in range(11):

    print(time)

    # send data to update
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    send_array[:,:,:,:] = 0
    send_array[1,:,porousStart:porousEnd,:] = Ub/eps
    # send_array[4,:,porousStart:porousEnd,:] = gradP
    CPL.send(send_array)

    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array)

# comm.Barrier()
CPL.finalize()
MPI.Finalize()