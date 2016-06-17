import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from cplpy import CPL
from draw_grid import draw_grid
from md_oo import MD

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
CFD_COMM = CPL.init(CPL.MD_REALM)
nprocs_realm = CFD_COMM.Get_size()

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
NProcs = np.product(npxyz)
xyzL = np.array([6.70820393, 17.88854382, 1.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

if (nprocs_realm != NProcs):
    print("Non-coherent number of processes in MD ", nprocs_realm,
            " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2])
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_md(cart_comm, xyzL, xyz_orig)

#Setup buffer to send CFD BC from MD
ncx = CPL.get("ncx"); dy = CPL.get("yl_cfd")/CPL.get("ncy")
ncy = np.floor(xyzL[1]/dy)
limits_CFD_BC = np.array([0, ncx, 0, 1, 0, 1], order='F', dtype=np.int32)
portion = CPL.my_proc_portion(limits_CFD_BC)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
A_send = np.zeros((2, ncxl, ncyl, nczl), order='F', dtype=np.float64)

#Setup buffer to recv constrained region
limits_MD_BC = np.array([0, ncx, 3, 4, 0, 1], order='F', dtype=np.int32)
portion = CPL.my_proc_portion(limits_MD_BC)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
A_recv = np.zeros((2, ncxl, ncyl, nczl), order='F', dtype=np.float64)

# Setup MD simulation object
md_cfd_dt_ratio = 25
dt = 0.005; Nsteps = 100000; tf = Nsteps*dt
time = np.arange(0.,tf,dt)
md = MD(dt=dt, wallwidth=[2.,0.], wallslide=[-1.,0.])

#Main run
for n,t in enumerate(time):

    print("MD time = ", md.tstep, md.time)

    # Calculate force
    md.force()

    #=======================================================
    # Call to CPL-LIBRARY
    # recieve u_CFD in constraint region
    # and force is applied
    # F = (1/tau)*(u_CFD - u_MD)
    #=======================================================
    if n%md_cfd_dt_ratio == 0:        
        A_recv, ierr = CPL.recv(A_recv, limits_MD_BC)
        u_CFD = A_recv[:,:,:,0]

    #Cell 7 is constrained
    md.constraint_force(u_CFD, 7)

    # Calculate velocity
    md.verlet()

    #Plot
    if n%md_cfd_dt_ratio == 0:        
        md.plot()

    #=======================================================
    #Call to CPL-LIBRARY to send u_MD at boundary
    #=======================================================
    if n%md_cfd_dt_ratio == 0:        
        u = md.get_velfield([ncx,ncy])
        #Cell 5 is sent
        A_send[0,:,0,0] = u[0,:,5]
        A_send[1,:,0,0] = u[1,:,5]
        CPL.send(A_send, limits_CFD_BC)




