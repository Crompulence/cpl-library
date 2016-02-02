from mpi4py import MPI
from cpl import CPL, get_olap_limits
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

#print "MD process " + str(realm_comm.Get_rank()) + " successfully initialised.\n"

# Scatter a 3 component array (MD processes are the receivers)
#scatter_array = np.zeros(3, order='F', dtype=np.float64)
"""
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
my_coords += 1# Fortran indices
extents = lib.proc_extents(my_coords, CPL.MD_REALM)
ncxl = extents[1] - extents[0] + 1  
ncyl = extents[3] - extents[2] + 1
nczl = extents[5] - extents[4] + 1
olap_limits = get_olap_limits(lib)
recv_array = -1*np.ones((3, ncxl, ncyl, nczl) , order='F', dtype=np.float64)
lib.scatter(scatter_array, olap_limits,recv_array)  

#if (np.all(recv_array > -1)):
#	print "MD rank: " + str(cart_comm.Get_rank()) + " recv_array: " + str(recv_array)

# Gather a 3 component array (MD processes are the senders)
gather_array = cart_rank * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
recv_array = np.zeros(3, order='F', dtype=np.float64)
lib.gather(gather_array, olap_limits, recv_array)
"""


