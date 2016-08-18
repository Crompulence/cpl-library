from __future__ import print_function, division
import sys
import cPickle

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

# Send dummy stress distribution (constant value of stress = 0) to MD
scatter_array = np.random.rand(9, cnstncx, cnstncy, cnstncz)


recv_array = np.zeros((9, 0, 0, 0), order='F', dtype=np.float64)
CPL.scatter(scatter_array, cnstFRegion, recv_array)

# Receive averaged velocities from LAMMPS socket
recv_array = np.zeros((4, velBCncx, velBCncy, velBCncz), order='F',
                      dtype=np.float64)
gather_array = np.zeros((4, 0, 0, 0), order='F', dtype=np.float64)
CPL.gather(gather_array, velBCRegion, recv_array)

dx, dy, dz = (Lx / NCx, Ly / NCy, Lz / NCz)
lines = ""
md_cell_coords = np.array(3)

dA = dx * dz
if (cnstFPortion[2] >= 0):
    for i in xrange(cnstFPortion[0], cnstFPortion[1]+1):
        for j in xrange(cnstFPortion[2], cnstFPortion[3]+1):
            for k in xrange(cnstFPortion[4], cnstFPortion[5] + 1):
                md_cell_coords = CPL.map_cell2coord(i, j, k)
                md_cell_coords = CPL.map_cfd2md_coord(md_cell_coords)
                [i_loc, j_loc, k_loc] = CPL.map_glob2loc_cell(cnstFPortion,
                                                              [i, j, k])
                lines += str(md_cell_coords[0] + dx/2.0) + " "\
                       + str(md_cell_coords[1] + dy/2.0) + " "\
                       + str(md_cell_coords[2] + dz/2.0) + " "\
                       + str(scatter_array[1, i_loc, j_loc, k_loc] * dA) + " "\
                       + str(scatter_array[4, i_loc, j_loc, k_loc] * dA) + " "\
                       + str(scatter_array[7, i_loc, j_loc, k_loc] * dA) + "\n"


lines = realm_comm.gather(lines, root=0)

myrank = realm_comm.Get_rank()
if myrank == 0:
    with open("cfd_forces.dat", "w") as forces_file:
        forces_file.writelines(lines)
