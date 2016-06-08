from mpi4py import MPI
from cplpy import CPL
import numpy as np

CPL = CPL()
nsteps = 1
dt = 0.2
density = 0.8

# Parameters of the cpu topology (cartesian grid)
NPx = 2
NPy = 2
NPz = 1
NProcs = NPx*NPy*NPz

# Parameters of the mesh topology (cartesian grid)
ncxyz = np.array([64, 18, 64], order='F', dtype=np.int32)
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "ERROR: Non-coherent number of processors."
    MPI.Abort(errorcode=1)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
CPL.setup_cfd(nsteps, dt, cart_comm, xyzL, xyz_orig, ncxyz, density)

cart_rank = cart_comm.Get_rank()
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
send_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)

for i in range(0, ncxl):
    for j in range(0, ncyl):
        for k in range(0, nczl):
            ii = i + portion[0]
            jj = j + portion[2]
            kk = k + portion[4]

            send_array[0, i, j, k] = ii
            send_array[1, i, j, k] = jj
            send_array[2, i, j, k] = kk

ierr = CPL.send(send_array, olap_limits)

MPI.COMM_WORLD.Barrier()

CPL.finalize()
MPI.Finalize()
