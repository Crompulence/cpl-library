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
xyzL = np.array([10.0, 16.0, 1.0], order='F', dtype=np.float64)
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
t0 = 0; tf = 30.; Nsteps = 10000
time = np.linspace(t0, tf, Nsteps)
dt = 0.005
md = MD(dt=dt, wallwidth=[2.,0.], wallslide=[-1.,0.])

#Main run
for n,t in enumerate(time):

    print("MD time = ", n,t)

    # Calculate force
    md.force()

    #=======================================================
    # Call to CPL-LIBRARY
    # recieve u_CFD in constraint region
    # and force is applied
    # F = (1/tau)*(u_CFD - u_MD)
    #=======================================================
    A_recv, ierr = CPL.recv(A_recv, limits_MD_BC)
    u_CFD = A_recv[0,:,:,:]
    print(A_recv.shape, u_CFD.shape)
    md.constraint_force(u_CFD, 8)

    # Calculate velocity
    md.verlet()

    #=======================================================
    #Call to CPL-LIBRARY to send u_MD at boundary
    #=======================================================
    u = md.get_velfield([ncx,ncy])
    A_send[0,:,0,0] = u[:,6,0] #Cell 6 is 
    A_send[1,:,0,0] = u[:,6,1]
    CPL.send(A_send, limits_CFD_BC)

    #Plot
    md.plot()
    


