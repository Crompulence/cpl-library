#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from CouetteAnalytical import CouetteAnalytical as CA

#Import CPL library
from cplpy import CPL

#initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([1., 1., 1.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

#Analytical solution
dt = 0.05
U = 1.
nu = 1e-02
Re = xyzL[1]*U/nu
ncx = CPL.get("ncx")
ncy = CPL.get("ncy")
ncz = CPL.get("ncz")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)
dx = xyzL[0]/float(ncx)
dy = xyzL[1]/float(ncy)
dz = xyzL[2]/float(ncz)
dV = dx*dy*dz
y = np.linspace(dy/2., xyzL[1] - dy/2., num=ncy)

#Main time loop
pcount = 0
for time in range(10000):

    print(time, recvbuf.shape)
    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf[...] = 0.
    CPL.send(sendbuf)

    ur = np.mean(recvbuf[0,:,:,:],(0,2))
    y_anal, u_anal = CAObj.get_vprofile(time*dt)

    #Plot values recieved from SediFOAM
    if time%100 == 0:
        plt.plot(ur, y,'go',ms=10,label="Recieved value from SediFOAM")
        plt.plot(u_anal, y_anal, 'k.-', label="Analytical Solution")
        plt.show()

CPL.finalize()
MPI.Finalize()




