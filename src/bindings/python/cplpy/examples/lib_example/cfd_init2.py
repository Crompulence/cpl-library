from mpi4py import MPI
from cplpy.cpl import CPL 
import numpy as np


#lib.set("output_mode", 0)
CPL = CPL()
nsteps = 1
dt = 0.2

# Parameters of the cpu topology (cartesian grid)
NPx = 3
NPy = 3
NPz = 3
NProcs = NPx*NPy*NPz
# Parameters of the mesh topology (cartesian grid)
ncxyz = np.array([18, 18, 18], order='F', dtype=np.int32)
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.create_comm(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "ERROR: Non-coherent number of processors."
    MPI.Abort(errorcode=1)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])

CPL.cfd_init2(nsteps, dt, cart_comm, xyzL, xyz_orig, ncxyz, 1.0)
ols = CPL.get_olap_limits()
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)

"""
#if (True):
#    print str(cart_comm.Get_rank())+" "+ "mycoords: " + str(my_coords)+ "extents:" + str(extents)+ " "+str(nsteps)+ str(dt)+ str(cart_comm)+ \
#    str(icoords)+ str(xyzL)+ str(ncxyz)+\
#    str(1.0)+ str(ijkcmax)+ str(ijkcmin)+ str(iTmin)+ str(iTmax)\
#    +str(jTmin)+ str(jTmax)+ str(kTmin)+ str(kTmax) 


[ncxl, ncyl, nczl] = ncxyz / [NPx, NPy, NPz]
#print "MD process " + str(realm_comm.Get_rank()) + " successfully initialised.\n"
# Scatter a 3 component array (CFD processes are the senders)
#print "cfdextents: " + str(extents) + "coord:" + str(my_coords)
#scatter_array = cart_rank * np.ones((3, ols[1]-ols[0]+1, ols[3]-ols[2]+1, ols[5]-ols[4]+1), order='F', dtype=np.float64)
cart_rank = cart_comm.Get_rank()
"""


"""
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
