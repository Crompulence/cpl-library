from __future__ import print_function
from mpi4py import MPI
from cplpy import CPL
import numpy as np


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
NPx = 4
NPy = 2
NPz = 2
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

cart_comm = MD_COMM.Create_cart([NPx, NPy, NPz])

CPL.setup_md(cart_comm, xyzL, xyz_orig)

# recv test
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)

recv_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
recv_array, ierr = CPL.recv(recv_array, olap_limits)

no_error = True
if CPL.overlap():
    rank = MD_COMM.Get_rank()
    for i in xrange(0, ncxl):
        for j in xrange(0, ncyl):
            for k in xrange(0, nczl):
                ii = i + portion[0]
                jj = j + portion[2]
                kk = k + portion[4]
                if (float(ii) - recv_array[0, i, j, k]) > 1e-8:
                    print("ERROR -- portion in x: %d %d " % (portion[0],
                          portion[1]) + " MD rank: %d " % rank +
                          " cell id: %d " % ii + " recv_array: %f" %
                          recv_array[0, i, j, k])
                    no_error = False

                if (float(jj) - recv_array[1, i, j, k]) > 1e-8:
                    print("ERROR -- portion in y: %d %d " % (portion[2],
                          portion[3]) + " MD rank: %d " % rank +
                          " cell id: %d " % jj + " recv_array: %f" %
                          recv_array[1, i, j, k])
                    no_error = False

                if (float(kk) - recv_array[2, i, j, k]) > 1e-8:
                    print("ERROR -- portion in z: %d %d " % (portion[4],
                          portion[5]) + " MD rank: %d " % rank +
                          " cell id: %d " % kk + " recv_array: %f" %
                          recv_array[2, i, j, k])
                    no_error = False


MD_COMM.Barrier()
if CPL.overlap() and no_error:
    print ("MD -- " + "(rank=%2d" % rank +
           ") CELLS HAVE BEEN RECEIVED CORRECTLY.\n", end="")
MPI.COMM_WORLD.Barrier()

#Free comms and finalise
MD_COMM.Free()
cart_comm.Free()

CPL.finalize()
MPI.Finalize()
