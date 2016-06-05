from __future__ import print_function, division
import cPickle
import sys

try:
    from mpi4py import MPI
    from cplpy import CPL
    import numpy as np
except ImportError as exc:
    print("ERROR: ", sys.exc_info()[0], exc, file=sys.stderr)
    MPI.COMM_WORLD.Abort(errorcode=1)


CPL = CPL()

CPL.set("output_mode", 1)


try:
    # Load parameters for the run
    params = cPickle.load(open("cfd_params.dic", "rb"))

    # Parameters of the cpu topology (cartesian grid)
    NPx = params["npx"]
    NPy = params["npy"]
    NPz = params["npz"]

    # Parameters of the mesh topology (cartesian grid)
    NCx = params["ncx"]
    NCy = params["ncy"]
    NCz = params["ncz"]

    # Parameters of the domain
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    MPI.COMM_WORLD.Abort(errorcode=1)

nsteps = 1
dt = 1.0
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
NProcs = NPx*NPy*NPz
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.COMM_WORLD.Abort(errorcode=1)

# Create cartesian communicator and initialize
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
CPL.setup_cfd(nsteps, dt, cart_comm, xyzL, xyz_orig, ncxyz, 1.0)

my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
olap_limits = CPL.get_olap_limits()

# Constrained region cell limits and number of cells
cnstFRegion = CPL.get_cnst_limits()
cnstFPortion = CPL.my_proc_portion(cnstFRegion)
[cnstncx, cnstncy, cnstncz] = CPL.get_no_cells(cnstFPortion)

# Velocity averaging region cell limits and number of cells
velBCRegion = np.copy(olap_limits)
velBCRegion[3] = velBCRegion[2]
velBCPortion = CPL.my_proc_portion(velBCRegion)
[velBCncx, velBCncy, velBCncz] = CPL.get_no_cells(velBCPortion)

# NOTE: Only 1 step for tests
for step in xrange(nsteps):
    # Send dummy stress distribution (constant value of stress = 0) to MD
    scatter_array = np.zeros((9, cnstncx, cnstncy, cnstncz), order='F',
                             dtype=np.float64)
    recv_array = np.zeros((9, 0, 0, 0), order='F', dtype=np.float64)
    CPL.scatter(scatter_array, cnstFRegion, recv_array)

    # Receive averaged velocities from LAMMPS socket
    recv_array = np.zeros((4, velBCncx, velBCncy, velBCncz), order='F',
                          dtype=np.float64)
    gather_array = np.zeros((4, 0, 0, 0), order='F', dtype=np.float64)
    CPL.gather(gather_array, velBCRegion, recv_array)

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
