import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from cplpy import CPL
from draw_grid import draw_grid
from cfd_oo import CFD

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
nprocs_realm = CFD_COMM.Get_size()

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
NProcs = np.product(npxyz)
xyzL = np.array([6.70820393, 17.88854382, 1.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([8, 8, 1], order='F', dtype=np.int32)

if (nprocs_realm != NProcs):
    print("Non-coherent number of processes in CFD ", nprocs_realm,
            " not equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2])
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz  )

#Setup buffer to get CFD BC from MD
ncx = CPL.get("ncx")
limits_CFD_BC = np.array([0, ncx, 0, 1, 0, 1], order='F', dtype=np.int32)
portion = CPL.my_proc_portion(limits_CFD_BC)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
A_recv = np.zeros((2, ncxl, ncyl, nczl), order='F', dtype=np.float64)

#Setup buffer to send constrained region
limits_MD_BC = np.array([0, ncx, 3, 4, 0, 1], order='F', dtype=np.int32)
portion = CPL.my_proc_portion(limits_MD_BC)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
A_send = np.zeros((2, ncxl, ncyl, nczl), order='F', dtype=np.float64)

#Set CFD simulation object
md_cfd_dt_ratio = 50
dt = 0.125; Nsteps = 100000/md_cfd_dt_ratio; tf = Nsteps*dt
time = np.arange(-dt,tf,dt)
uwall = 1.
cfd = CFD(nu=0.575, dt=dt, 
          xsize = ncxyz[0], ysize = ncxyz[1]+2,
          xmin = xyz_orig[0], xmax = xyzL[0],
          ymin = xyz_orig[1], ymax = xyzL[1])

#Main Run
for n,t in enumerate(time):

    print("CFD time = ", n,t)

    #===============================================
    # Call to CPL-LIBRARY goes here to
    # send u_CFD in constraint region
    #===============================================
    A_send[0,:,0,0] = cfd.u[:,2]
    CPL.send(A_send, limits_MD_BC)

    #===============================================
    # Call to CPL-LIBRARY goes here to
    # recieve u_MD to set bottom boundary
    #===============================================
    umd, ierr = CPL.recv(A_recv, limits_CFD_BC)
    bottomwall = umd[0,:,0,0]
    cfd.set_bc(topwall=uwall, bottomwall=bottomwall)

    #plot
    cfd.plot()

    #Update CFD timestep
    cfd.update_time()


