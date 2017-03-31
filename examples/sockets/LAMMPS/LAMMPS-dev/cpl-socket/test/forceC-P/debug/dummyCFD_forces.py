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


cpllib = CPL()

cpllib.set("output_mode", 1)

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
realm_comm = cpllib.init(cpllib.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processes is not coherent.", file=sys.stderr)
    MPI.COMM_WORLD.Abort(errorcode=1)

# Create cartesian communicator and initialize
cpllib.set_timing(0, nsteps, dt)
cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cpllib.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

my_coords = cart_comm.Get_coords(cart_comm.Get_rank())
my_coords = np.array(my_coords, order='F', dtype=np.int32)
olap_limits = cpllib.get_olap_limits()

# Constrained region cell limits and number of cells
cnstFRegion = cpllib.get_cnst_limits()
cnstFPortion = cpllib.my_proc_portion(cnstFRegion)
[cnstncx, cnstncy, cnstncz] = cpllib.get_no_cells(cnstFPortion)

# Velocity averaging region cell limits and number of cells
velBCRegion = np.copy(olap_limits)
velBCRegion[3] = velBCRegion[2]
velBCPortion = cpllib.my_proc_portion(velBCRegion)
[velBCncx, velBCncy, velBCncz] = cpllib.get_no_cells(velBCPortion)

# Send dummy random stress distribution (constant value of stress = 0) to MD
send_array = 10*np.ones((9, cnstncx, cnstncy, cnstncz), order='F', dtype=np.float64)

cpllib.send(send_array, cnstFRegion)

# Receive averaged velocities from LAMMPS socket
recv_array = np.zeros((4, velBCncx, velBCncy, velBCncz), order='F',
                      dtype=np.float64)
cpllib.recv(recv_array, velBCRegion)

dx, dy, dz = (Lx / NCx, Ly / NCy, Lz / NCz)
dA = dx * dz

# Dump 1, 4 an 7 components of the stress tensor multiplied by the area of the cell
func = lambda x: x*dA
cpllib.dump_region(cnstFRegion, send_array, "cfd_forces.dat", realm_comm,                                                                                                                                                            
        components={1:func, 4:func, 7:func}, coords="other")

cpllib.finalize()
