from __future__ import print_function, division
from mpi4py import MPI
from cplpy.cpl import CPL
import numpy as np
import cPickle
import sys

comm_world = MPI.COMM_WORLD
CPL = CPL()

# Do not show any info to the stdin
CPL.set("output_mode", 0)

nsteps = 1
dt = 0.2

# Load parameters for the run
params = cPickle.load(open("cfd_params.dic", "rb"))

# Parameters of the cpu topology (cartesian grid)
try:
    NPx = params["npx"]
    NPy = params["npy"]
    NPz = params["npz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)

# Parameters of the mesh topology (cartesian grid)
try:
    NCx = params["ncx"]
    NCy = params["ncy"]
    NCz = params["ncz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)

# Parameters of the domain
try:
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)

NProcs = NPx*NPy*NPz
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)


# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    comm_world.Abort(errorcode=1)

# Create cartesian communicator and build lookup table
# for the grid topology.
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])

CPL.setup_cfd(nsteps, dt, cart_comm, xyzL, xyz_orig, ncxyz, 1.0)


# Sending cell coordinates from CFD to MD

olap_region = CPL.get_olap_region_region()
portion = CPL.my_proc_portion(olap_region)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
send_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
for iglob in xrange(portion[0], portion[1]):
    for jglob in xrange(portion[2], portion[3]):
        for kglob in xrange(portion[4], portion[5]):
            iloc, jloc, kloc, = CPL.map_glob2loc(portion, [iglob, jglob, kglob])
            send_array[0:3, iloc, jloc, kloc] = [iglob, jglob, kglob]

recv_array = np.zeros((3, 0, 0, 0), order='F', dtype=np.float64)
CPL.scatter(send_array, olap_region, recv_array)


recv_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
send_array = np.zeros((3, 0, 0, 0), order='F', dtype=np.float64)
CPL.gather(send_array, olap_region, recv_array)
for icfd in xrange(portion[0], portion[1]):
    for jcfd in xrange(portion[2], portion[3]):
        for kcfd in xrange(portion[4], portion[5]):
            iloc, jloc, kloc, = CPL.map_glob2loc(portion, [icfd, jcfd, kcfd])
            [imd, jmd, kmd] = recv_array[0:3, iloc, jloc, kloc]
            lines += str(imd) + str(jmd) + str(kmd) + str(jmd


lines = realm_comm.gather(lines, root=0)
myrank = realm_comm.Get_rank()
if myrank == 0:
  with open("cfd_recv_cells.dat", "w") as cells_file:
      cells_file.writelines(lines)
