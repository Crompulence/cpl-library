from mpi4py import MPI
from cplpy.cpl import CPL 
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
limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
send_array = np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)

for i in range(portion[0],portion[1]+1):
    for j in range(portion[2],portion[3]+1):
        for k in range(portion[4],portion[5]+1):
            ii = i - portion[0]
            jj = j - portion[2]
            kk = k - portion[4]

            send_array[0,ii,jj,kk] = i+1
            send_array[1,ii,jj,kk] = j+1
            send_array[2,ii,jj,kk] = k+1

ierr = CPL.send(send_array, limits)

MPI.COMM_WORLD.Barrier()
