from __future__ import print_function, division
from mpi4py import MPI
from cplpy.cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits
import numpy as np
import cPickle
import sys


lib = CPL()

lib.set("output_mode", 0)

nsteps = 100
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
    MPI.Abort()

# Parameters of the mesh topology (cartesian grid)
try:
    NCx = params["ncx"]
    NCy = params["ncy"]
    NCz = params["ncz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    MPI.Abort()

# Parameters of the domain
try:
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    MPI.Abort()

NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
dx = Lx / NCx
dy = Ly / NCy
dz = Lz / NCz
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)


# Create communicators and check that number of processors is consistent
realm_comm = lib.create_comm(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.Abort()

# Create mesh grid points (CPL format)
xg, yg, zg = create_CPL_cart_3Dgrid(NCx, NCy, NCz, dx, dy, dz)

# Partition of the number of processors
ijkcmax = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
ijkcmin = np.array([1, 1, 1], order='F', dtype=np.int32)


# Create cartesian communicator and build lookup table
# for the grid topology.
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cart_rank = cart_comm.Get_rank()
icoords = np.zeros((3, NProcs), order='F', dtype=np.int32)


iTmin = np.zeros(NPx, order='F', dtype=np.int32)
iTmax = np.zeros(NPx, order='F', dtype=np.int32)
jTmin = np.zeros(NPy, order='F', dtype=np.int32)
jTmax = np.zeros(NPy, order='F', dtype=np.int32)
kTmin = np.zeros(NPz, order='F', dtype=np.int32)
kTmax = np.zeros(NPz, order='F', dtype=np.int32)


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

lib.cfd_init(nsteps, dt, cart_comm, icoords, npxyz, xyzL, ncxyz,
             1.0, ijkcmax, ijkcmin, iTmin, iTmax, jTmin, jTmax, kTmin,
             kTmax, xg, yg, zg)

olap_limits = get_olap_limits(lib)
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
# Fortran indices
my_coords += 1
extents = lib.proc_extents(my_coords, CPL.CFD_REALM)
