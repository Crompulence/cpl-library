from __future__ import print_function, division
from mpi4py import MPI
from cplpy.cpl import CPL, get_olap_limits
import numpy as np
import cPickle
import sys

comm = MPI.COMM_WORLD
lib = CPL()

# Do not show any info to the stdin
lib.set("output_mode", 0)

dt = 0.1

# Load parameters for the run
params = cPickle.load(open("md_params.dic", "rb"))

# Parameters of the cpu topology (cartesian grid)
try:
    NPx = params["npx"]
    NPy = params["npy"]
    NPz = params["npz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    MPI.Abort()

NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Parameters of the domain
try:
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    MPI.Abort()

global_domain = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = lib.create_comm(CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.Abort()

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cart_rank = cart_comm.Get_rank()
icoords = np.zeros((3, NProcs), order='F', dtype=np.int32)

for rank in xrange(NProcs):
    cart_coords = cart_comm.Get_coords(rank)
    icoords[:, rank] = cart_coords
icoords += 1
nsteps, initialstep = \
        lib.md_init(dt, cart_comm, icoords, npxyz, global_domain, 1)

# Scatter a 3 component array (MD processes are the receivers)
my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
# Fortran indices
my_coords += 1
olap_limits = get_olap_limits(lib)
# extents = lib.proc_portion(my_coords, CPL.MD_REALM, olap_limits)
extents = lib.proc_extents(my_coords, CPL.MD_REALM)
ncxl = extents[1] - extents[0] + 1
ncyl = extents[3] - extents[2] + 1
nczl = extents[5] - extents[4] + 1
