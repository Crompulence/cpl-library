from __future__ import print_function, division
from mpi4py import MPI
from cplpy.cpl import CPL
import numpy as np
import cPickle
import sys


lib = CPL()

lib.set("output_mode", 1)


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

nsteps = 1
dt = 1.0
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
NProcs = NPx*NPy*NPz
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = lib.create_comm(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.Abort()

# Create cartesian communicator and initialize
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
lib.cfd_init2(nsteps, dt, cart_comm, xyzL, xyz_orig, ncxyz, 1.0)

my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
olap_limits = lib.get_olap_limits()

# Constrained region cell limits and number of cells
cnstFRegion = lib.get_cnst_limits()
cnstFPortion = lib.my_proc_portion(cnstFRegion)
[cnstncx, cnstncy, cnstncz] = lib.get_no_cells(cnstFPortion)

# Velocity averaging region cell limits and number of cells
velBCRegion = np.copy(olap_limits)
velBCRegion[3] = velBCRegion[2]
velBCPortion = lib.my_proc_portion(velBCRegion)
[velBCncx, velBCncy, velBCncz] = lib.get_no_cells(velBCPortion)

# Send dummy stress distribution (constant value of stress = 0) to MD
scatter_array = np.zeros((9, cnstncx, cnstncy, cnstncz), order='F',
                         dtype=np.float64)
recv_array = np.zeros((9, 0, 0, 0), order='F', dtype=np.float64)
lib.scatter(scatter_array, cnstFRegion, recv_array)

# Receive averaged velocities from LAMMPS socket
recv_array = np.zeros((4, velBCncx, velBCncy, velBCncz), order='F', dtype=np.float64)
gather_array = np.zeros((4, 0, 0, 0), order='F', dtype=np.float64)
lib.gather(gather_array, velBCRegion, recv_array)
dx, dy, dz = (Lx / NCx, Ly / NCy, Lz / NCz)
lines = ""

if (velBCPortion[2] >= 0):
    for i in xrange(velBCncx):
        for k in xrange(velBCncz):
            icoord = my_coords[0] * velBCncx + i
            kcoord = my_coords[2] * velBCncz + k
            lines += str(icoord*dx + dx/2.0) + " "\
                   + str(kcoord*dz + dz/2.0) + " "\
                   + str(recv_array[0, i, 0, k]) + " "\
                   + str(recv_array[1, i, 0, k]) + " "\
                   + str(recv_array[2, i, 0, k]) + "\n"

lines = realm_comm.gather(lines, root=0)

myrank = realm_comm.Get_rank()
if myrank == 0:
    with open("cfd_vels.dat", "w") as vels_file:
        vels_file.writelines(lines)
