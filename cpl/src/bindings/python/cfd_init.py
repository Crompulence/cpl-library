from mpi4py import MPI
from cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits
import numpy as np



lib = CPL()

nsteps = 100
dt = 0.2

# Parameters of the cpu topology (cartesian grid)
NPx = 3
NPy = 3
NPz = 3
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)
# Parameters of the mesh topology (cartesian grid)
NCx = 18
NCy = 18
NCz = 18
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
Lx = 10.0
Ly = 10.0
Lz = 10.0
dx = Lx / NCx
dy = Ly / NCy
dz = Lz / NCz
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)


# Create communicators and check that number of processors is consistent
realm_comm = lib.create_comm(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "Non-coherent number of processes"
    exit()
# Create mesh grid points (CPL format)
xg, yg, zg = create_CPL_cart_3Dgrid(NCx, NCy, NCz, dx, dy , dz)

# Partition of the number of processors
ijkcmax = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
ijkcmin = np.array([1, 1, 1], order='F', dtype=np.int32)


#Minimun and maximun cell number for each process
"""
[ncxl, ncyl, nczl] = ncxyz / npxyz
iTmin = np.arange(1, ncxyz[0], ncxl, dtype=np.int32)
iTmax = iTmin + ncxl - 1
jTmin = np.arange(1, ncxyz[1], ncyl, dtype=np.int32)
jTmax = jTmin + ncyl - 1
kTmin = np.arange(1, ncxyz[2], nczl, dtype=np.int32)
kTmax = kTmin + nczl - 1
"""

# Create cartesian communicator and build lookup table
# for the grid topology.
#cart_comm = cart_create(realm_comm, [NCx, NCy, NCz], [0,0,0], 
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cart_rank = cart_comm.Get_rank()
icoords = np.zeros((3, NProcs), order='F', dtype=np.int32)


iTmin = np.zeros( NProcs, order='F', dtype=np.int32)
iTmax = np.zeros( NProcs, order='F', dtype=np.int32)
jTmin = np.zeros( NProcs, order='F', dtype=np.int32)
jTmax = np.zeros( NProcs, order='F', dtype=np.int32)
kTmin = np.zeros( NProcs, order='F', dtype=np.int32)
kTmax = np.zeros( NProcs, order='F', dtype=np.int32)


	
[ncxl, ncyl, nczl] = ncxyz / npxyz
for rank in xrange(NProcs):
	cart_coords = cart_comm.Get_coords(rank)
	icoords[:, rank] = cart_coords
	x = cart_coords[0]
	y = cart_coords[1]
	z = cart_coords[2]
	iTmin[rank] = x*ncxl + 1
	iTmax[rank] = iTmin[rank] + ncxl - 1
	jTmin[rank] = y*ncyl + 1
	jTmax[rank] = jTmin[rank] + ncyl - 1
	kTmin[rank] = z*nczl + 1
	kTmax[rank] = kTmin[rank] + nczl - 1

# Fortran indexing
icoords += 1
if cart_comm.Get_rank() == 0:
    print "iTmin: " + str(iTmin)
    print "iTmax: " + str(iTmax)
    print "jTmin: " + str(jTmin)
    print "jTmax: " + str(jTmax)
    print "kTmin: " + str(kTmin)
    print "kTmax: " + str(kTmax)
    print "ijkcmin: " + str(ijkcmin)
    print "ijkcmax: " + str(ijkcmax)
    print "topo_coords: " + str(icoords)
    print "xg: " +str(xg) + "len: " + str(len(xg))
    print "yg: " +str(yg) + "len: " + str(len(yg))
    print "zg: " +str(zg) + "len: " + str(len(zg))



lib.cfd_init(nsteps, dt, cart_comm, icoords, npxyz, xyzL, ncxyz,
            1.0, ijkcmax, ijkcmin, iTmin, iTmax, jTmin, jTmax, kTmin,
            kTmax, xg, yg, zg)


#print "MD process " + str(realm_comm.Get_rank()) + " successfully initialised.\n"

"""
# Scatter a 3 component array (CFD processes are the senders)
scatter_array = cart_rank * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
olap_limits = get_olap_limits(lib)
recv_array = np.zeros(3 , order='F', dtype=np.float64)
#lib.scatter(scatter_array, olap_limits, recv_array)  

if (np.all(recv_array > -1)):
	print "CFD rank: " + str(cart_comm.Get_rank())# + " recv_array: " + str(recv_array)
# Gather a 3 component array (CFD processes are the receivers)
recv_array = -1 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
gather_array = np.zeros(0, order='F', dtype=np.float64)
lib.gather(gather_array, olap_limits, recv_array)


print "something"
"""



#if (comm.Get_rank() == 1):
#    print "kTmin: ", kTmin, " kTmax: ",  kTmax
#    print xg
