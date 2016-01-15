from mpi4py import MPI
from cpl import CPL
import numpy as np
from ctypes import c_int

comm = MPI.COMM_WORLD
lib = CPL()

nsteps = c_int(100)
initialstep = c_int(0)
dt = 0.1

# Parameters of the cpu topology (cartesian grid)
NPx = 3
NPy = 3
NPz = 3
NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Domain topology
Lx = 10.0
Ly = 10.0
Lz = 10.0
global_domain = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)

realm_comm = lib.create_comm(CPL.MD_REALM)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])
cart_rank = cart_comm.Get_rank()
icoords = np.zeros((NProcs, 3), order='F', dtype=np.int32)

lib.md_init(nsteps, initialstep, dt, cart_comm, icoords, npxyz, global_domain, 1.0)




