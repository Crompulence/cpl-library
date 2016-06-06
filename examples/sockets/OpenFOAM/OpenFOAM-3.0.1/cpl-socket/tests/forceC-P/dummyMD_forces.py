from __future__ import print_function, division
from mpi4py import MPI
from cplpy import CPL
import numpy as np
import cPickle
import sys


lib = CPL()

lib.set("output_mode", 1)


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

dt = 1.0
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
NProcs = NPx*NPy*NPz
ncxyz = np.array([NCx, NCy, NCz], order='F', dtype=np.int32)
xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = lib.init (CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.Abort()

# Create cartesian communicator and initialize
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
nsteps, initial_step = lib.setup_md (dt, cart_comm, xyzL, xyz_orig, 1.0)
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

for step in xrange(nsteps):
    # Send dummy stress distribution (constant value of stress = 0) to MD
    send_array = np.zeros((9, 0, 0, 0), order='F', dtype=np.float64)
    recv_array = np.ones((9, cnstncx, cnstncy, cnstncz), order='F', dtype=np.float64)
    lib.scatter(send_array, cnstFRegion, recv_array)
    #if (cnstFPortion[2] >= 0):
    #    print (recv_array)
#with open("forces.dat", "w") as forces_file:
#    forces_file.writelines( str (recv_array))
# Receive averaged velocities from LAMMPS socket
    #send_array = 0.5*np.ones((4, velBCncx, velBCncy, velBCncz), order='F', dtype=np.float64)
    send_array = np.random.rand(4, velBCncx, velBCncy, velBCncz)
    recv_array = np.zeros((4, 0, 0, 0), order='F', dtype=np.float64)
    lib.gather(send_array, velBCRegion, recv_array)
    """
dx, dy, dz = (Lx / NCx, Ly / NCy, Lz / NCz)
lines = ""
md_cell_coords = np.array(3)

dA = dx * dz
if (cnstFPortion[2] >= 0): 
    for i in xrange(cnstFPortion[0], cnstFPortion[1] + 1):
        for j in xrange(cnstFPortion[2], cnstFPortion[3] + 1 ):
            for k in xrange(cnstFPortion[4], cnstFPortion[5] + 1):
                md_cell_coords = lib.grid_cell2xyz(i, j, k) 
                md_cell_coords = lib.map_cfd2md_global(md_cell_coords)
                [i_loc, j_loc, k_loc] = lib.cell_glob2loc(cnstFPortion, [i, j, k])
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
"""
