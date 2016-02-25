from mpi4py import MPI
from cplpy.cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits
import numpy as np



lib = CPL()

#lib.set("output_mode", 0)

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
#print "xg:" + str(xg) + "yg:" + str(yg) + "zg:" + str(zg)

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


iTmin = np.zeros( NPx, order='F', dtype=np.int32)
iTmax = np.zeros( NPx, order='F', dtype=np.int32)
jTmin = np.zeros( NPy, order='F', dtype=np.int32)
jTmax = np.zeros( NPy, order='F', dtype=np.int32)
kTmin = np.zeros( NPz, order='F', dtype=np.int32)
kTmax = np.zeros( NPz, order='F', dtype=np.int32)


	
[ncxl, ncyl, nczl] = ncxyz / npxyz
for rank in xrange(NProcs):
	cart_coords = cart_comm.Get_coords(rank)
	icoords[:, rank] = cart_coords
	x = cart_coords[0]
	y = cart_coords[1]
	z = cart_coords[2]
	iTmin[x] = x*ncxl + 1
	iTmax[x] = iTmin[x] + ncxl - 1
	jTmin[y] = y*ncyl + 1
	jTmax[y] = jTmin[y] + ncyl - 1
	kTmin[z] = z*nczl + 1
	kTmax[z] = kTmin[z] + nczl - 1


# Fortran indexing
icoords += 1
"""
print "iTmin: " + str(iTmin)
print "iTmax: " + str(iTmax)
print "jTmin: " + str(jTmin)
print "jTmax: " + str(jTmax)
print "kTmin: " + str(kTmin)
print "kTmax: " + str(kTmax)
print "ijkcmin: " + str(ijkcmin)
print "ijkcmax: " + str(ijkcmax)
"""
"""
print "topo_coords: " + str(icoords)
print "xg: " +str(xg) + "len: " + str(len(xg))
print "yg: " +str(yg) + "len: " + str(len(yg))
print "zg: " +str(zg) + "len: " + str(len(zg))
"""


lib.cfd_init(nsteps, dt, cart_comm, icoords, npxyz, xyzL, ncxyz,
            1.0, ijkcmax, ijkcmin, iTmin, iTmax, jTmin, jTmax, kTmin,
            kTmax, xg, yg, zg)
#if (cart_comm.Get_rank() == 2):
ols = get_olap_limits(lib)
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
my_coords += 1# Fortran indices
extents = lib.proc_extents(my_coords, CPL.CFD_REALM)

#if (True):
#    print str(cart_comm.Get_rank())+" "+ "mycoords: " + str(my_coords)+ "extents:" + str(extents)+ " "+str(nsteps)+ str(dt)+ str(cart_comm)+ \
#    str(icoords)+ str(npxyz)+ str(xyzL)+ str(ncxyz)+\
#    str(1.0)+ str(ijkcmax)+ str(ijkcmin)+ str(iTmin)+ str(iTmax)\
#    +str(jTmin)+ str(jTmax)+ str(kTmin)+ str(kTmax) 


#print "MD process " + str(realm_comm.Get_rank()) + " successfully initialised.\n"
# Scatter a 3 component array (CFD processes are the senders)
"""
#print "cfdextents: " + str(extents) + "coord:" + str(my_coords)
#scatter_array = cart_rank * np.ones((3, ols[1]-ols[0]+1, ols[3]-ols[2]+1, ols[5]-ols[4]+1), order='F', dtype=np.float64)
scatter_array = cart_rank * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
#print ols , "len: ", scatter_array.size
recv_array = np.zeros(3 , order='F', dtype=np.float64)
lib.scatter(scatter_array, ols, recv_array)  
"""
"""
# Gathe a 3 component array (CFD processes are the receivers)
recv_array = -1 * np.ones((3, ols[1], ols[3], ols[5]), order='F', dtype=np.float64)
#recv_array = -1 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
gather_array = np.zeros((3,0,0,0), order='F', dtype=np.float64)
lib.gather(gather_array, ols, recv_array)

if (np.any(recv_array > -1)):
    print "CFD rank: " + str(cart_comm.Get_rank())# + " recv_array: " + str(recv_array)

"""
#if (comm.Get_rank() == 1):
#    print "kTmin: ", kTmin, " kTmax: ",  kTmax
#    print xg
# Send test
"""
ols = get_olap_limits(lib)
send_array = cart_rank * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
lib.send(send_array, ols[0], ols[1], ols[2], ols[3], ols[4], ols[5])
"""
"""
ols = get_olap_limits(lib)
recv_array = -1 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
lib.recv(recv_array, ols[0], ols[1], ols[2], ols[3], ols[4], ols[5])
if (np.all(recv_array > -1)):
    print "mycoords: " + str(my_coords)+ "CFD rank: " + str(cart_comm.Get_rank()) + " recv_array: " + str(recv_array)
"""
