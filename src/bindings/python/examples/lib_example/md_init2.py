from mpi4py import MPI
from cplpy.cpl import CPL
import numpy as np

comm = MPI.COMM_WORLD
CPL = CPL()

#lib.set("output_mode", 0)

dt = 0.1

# Parameters of the cpu topology (cartesian grid)
NPx = 3
NPy = 3
NPz = 3
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Domain topology
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.create_comm(CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "Non-coherent number of processes"
    MPI.Abort(errorcode=1)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])

nsteps, initialstep = CPL.md_init2(dt, cart_comm, xyzL, xyz_orig, 1.0)
"""
# Scatter a 3 component array (MD processes are the receivers)
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
#@my_coords += 1# Fortran indices
olap_limits = get_olap_limits(lib)
#extents = lib.proc_portion(my_coords, CPL.MD_REALM, olap_limits)
extents = lib.proc_extents(my_coords, CPL.MD_REALM)
#print "mdextents: " + str(extents) + "coord:" + str(my_coords)
ncxl = extents[1] - extents[0] + 1  
ncyl = extents[3] - extents[2] + 1
nczl = extents[5] - extents[4] + 1
"""
# if (True):
#     print str(cart_comm.Get_rank())+" "+ "mycoords: " + str(my_coords)+ "extents:" + str(extents)+ " "+str(dt)+ str(cart_comm)+ \
#     str(icoords)+ str(npxyz)+ str(global_domain)
#

"""
scatter_array = np.zeros(3, order='F', dtype=np.float64)
recv_array = -1*np.ones((3, ncxl, ncyl, nczl) , order='F', dtype=np.float64)
lib.scatter(scatter_array, olap_limits,recv_array)  

if (np.any(recv_array > -1)):
    print "extents:" + str(extents) + "mycoords: "+ str(my_coords) + "MD rank: " + str(cart_comm.Get_rank()) + " recv_array: " + str(recv_array)
"""
"""
# Gather a 3 component array (MD processes are the senders)
gather_array = 5 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
recv_array = np.zeros(3, order='F', dtype=np.float64)
lib.gather(gather_array, olap_limits, recv_array)
"""
"""
# Send test
ols = get_olap_limits(lib)
recv_array = -1 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
lib.recv(recv_array, ols[0], ols[1], ols[2], ols[3], ols[4], ols[5])

if (np.any(recv_array > -1)):
	print "extents:" + str(extents) + "mycoords: "+ str(my_coords) + "MD rank: " + str(cart_comm.Get_rank()) + " recv_array: " + str(recv_array)
"""
"""
ols = get_olap_limits(lib)
send_array = cart_rank * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
lib.send(send_array, ols[0], ols[1], ols[2], ols[3], ols[4], ols[5])
"""
