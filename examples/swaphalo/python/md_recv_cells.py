from __future__ import print_function
from mpi4py import MPI
from cplpy import CPL
import numpy as np
import matplotlib.pyplot as plt

def read_input(filename):
    with open(filename, 'r') as f:
        content = f.read()

    dic = {}
    for i in content.split('\n'):
        if i.find("!") != -1:
            name = i.split("!")[1]
            value = i.split('!')[0].replace(' ', '')
            dic[name] = value
    return dic

comm = MPI.COMM_WORLD
CPL = CPL()

# Parameters of the cpu topology (cartesian grid)
dt = 0.1
NPx = 1
NPy = 1
NPz = 1
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Domain topology
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
MD_COMM = CPL.init(CPL.MD_REALM)
nprocs_realm = MD_COMM.Get_size()

if (nprocs_realm != NProcs):
    print("Non-coherent number of processes")
    MPI.Abort(errorcode=1)

cart_comm = MD_COMM.Create_cart([NPx, NPy, NPz], periods=[True,True,True])

CPL.setup_md(cart_comm, xyzL, xyz_orig)

# recv test
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)

nh = [1,1,1]
A = np.zeros((3, ncxl+2*nh[0], ncyl+2*nh[1], nczl+2*nh[2]), order='F', dtype=np.float64)

for i in range(nh[0], ncxl+nh[0]):
    for j in range(nh[1], ncyl+nh[1]):
        for k in range(nh[2], nczl+nh[2]):
            ii = i + portion[0]
            jj = j + portion[2]
            kk = k + portion[4]

            A[0, i, j, k] = ii + jj*ncxl*NPx
            A[1, i, j, k] = 0.4
            A[2, i, j, k] = 0.5

plt.imshow(A[0,:,:,3], interpolation="None")
plt.colorbar()
plt.show()

A = CPL.swaphalos(A)

plt.imshow(A[0,:,:,3], interpolation="None")
plt.colorbar()
plt.show()


#Free comms and finalise
MD_COMM.Free()
cart_comm.Free()

CPL.finalize()
MPI.Finalize()
