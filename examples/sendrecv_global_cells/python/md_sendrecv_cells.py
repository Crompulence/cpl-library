from mpi4py import MPI
from cplpy.cpl import CPL
import numpy as np

def read_input(filename):
    with open(filename,'r') as f:
        content = f.read()

    dic = {}
    for i in content.split('\n'):
        if i.find("!") != -1:
            name = i.split("!")[1]
            value = i.split('!')[0].replace(' ','')
            dic[name] = value
    
    return dic

comm = MPI.COMM_WORLD
CPL = CPL()

# Parameters of the cpu topology (cartesian grid)
dt = 0.1
NPx = 4
NPy = 2
NPz = 2
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Domain topology
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print "Non-coherent number of processes"
    MPI.Abort(errorcode=1)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])

nsteps, initialstep = CPL.setup_md(dt, cart_comm, xyzL, xyz_orig, 1.0)

# Send test
ols = CPL.get_olap_limits()
portion = CPL.my_proc_portion(ols)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)

recv_array = -1 * np.ones((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
recv_array, ierr = CPL.recv(recv_array, ols)

if (np.any(recv_array > -1)):
	print "extents:" + str(portion) + "MD rank: " + str(cart_comm.Get_rank()) + " recv_array: " + str(recv_array)

