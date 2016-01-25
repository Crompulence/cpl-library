from mpi4py import MPI
from cpl import CPL
import numpy as np

comm = MPI.COMM_WORLD
lib = CPL()

dt = 0.1

# Parameters of the cpu topology (cartesian grid)
NPx = 3
NPy = 3
NPz = 3
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Domain topology
Lx = 10.0
Ly = 10.0
Lz = 10.0
global_domain = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = lib.create_comm(CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "Non-coherent number of processes"
    exit()

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cart_rank = cart_comm.Get_rank()
icoords = np.zeros((3, NProcs), order='F', dtype=np.int32)

for rank in xrange(NProcs):
    cart_coords = cart_comm.Get_coords(rank)
    icoords[:, rank] = cart_coords 
icoords += 1
nsteps, initialstep = lib.md_init(dt, cart_comm, icoords, npxyz, global_domain, 1)
scatter_array = np.
lib.scatter

print "MD process " + str(realm_comm.Get_rank()) + " successfully initialised.\n"
