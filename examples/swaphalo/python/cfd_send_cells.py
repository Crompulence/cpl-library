from mpi4py import MPI
from cplpy import CPL
import numpy as np


CPL = CPL()
nsteps = 1
dt = 0.2
density = 0.8

# Parameters of the cpu topology (cartesian grid)
NPx = 1
NPy = 1
NPz = 1
NProcs = NPx*NPy*NPz

# Parameters of the mesh topology (cartesian grid)
ncxyz = np.array([8, 8, 8], order='F', dtype=np.int32)
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
CFD_COMM = CPL.init(CPL.CFD_REALM)
nprocs_realm = CFD_COMM.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: Non-coherent number of processors.")
    MPI.Abort(errorcode=1)

cart_comm = CFD_COMM.Create_cart([NPx, NPy, NPz])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

cart_rank = cart_comm.Get_rank()
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)

nh = [1,1,1]
A = np.zeros((3, ncxl+2*nh[0], ncyl+2*nh[1], nczl+2*nh[2]), order='F', dtype=np.float64)

for i in range(0, ncxl):
    for j in range(0, ncyl):
        for k in range(0, nczl):
            ii = i + portion[0]+nh[0]
            jj = j + portion[2]+nh[1]
            kk = k + portion[4]+nh[2]

            A[0, ii, jj, kk] = 0.5
            A[1, ii, jj, kk] = 0.5
            A[2, ii, jj, kk] = 0.5

#A = CPL.swaphalos(A)


MPI.COMM_WORLD.Barrier()

CFD_COMM.Free()
cart_comm.Free()

CPL.finalize()
MPI.Finalize()
