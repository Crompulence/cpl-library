#!/usr/bin/env python

import numpy as np
from mpi4py import MPI
from CouetteAnalytical import CouetteAnalytical as CA
import matplotlib.pyplot as plt

# From /cpl-library/utils/
from plot_parallel import collect_data

#Import CPL library
from cplpy import CPL

plot_stuff = False

#initialise MPI
comm = MPI.COMM_WORLD
root = 0

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([1., 1., 1.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
cart_COMM = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_md(cart_COMM, xyzL, xyz_orig)
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

#Analytical solution
dt = 0.05
U = 1.
nu = 1.004e-2
Re = xyzL[1]*U/nu
ncx = CPL.get("ncx")
ncy = CPL.get("ncy")
ncz = CPL.get("ncz")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=200*ncy)
dx = xyzL[0]/float(ncx)
dy = xyzL[1]/float(ncy)
dz = xyzL[2]/float(ncz)
dV = dx*dy*dz
y = np.linspace(dy/2., xyzL[1] - dy/2., num=ncy)

#Main time loop
pcount = 0
error_list = []
for time in range(2000):

    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    #Collect data for plotting while next loop evolves in CFD
    globalrecvbuf = collect_data(recvbuf, cart_COMM, [ncx, ncy, ncz], plotrank=root)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf[...] = 0.
    CPL.send(sendbuf)

    if MD_COMM.rank == root:

        ur = np.mean(globalrecvbuf[0,:,:,:],(0,2))
        t = max((time-.5)*dt,0.)
        y_anal, u_anal = CAObj.get_vprofile(t)
        ua = u_anal[1:-1:2]
        error = (ua - ur)/U
        if time%20 == 0:
            print("Time= " + str(time) + " Error_vs_analytical_soln.= "
                           + str(100.*np.sqrt(np.sum(error**2))) + " %")

        error_list.append(100.*np.sqrt(np.sum(error**2)))

        if plot_stuff:
            
            #Plot values recieved from SediFOAM
            if time%10000 == 0:
                plt.plot(ur, y,'go',ms=10,label="Recieved value from SediFOAM")
                plt.plot(u_anal, y_anal, 'k.-', label="Analytical Solution")

                #Plot Error
                plt.plot(100.*error, y)
                plt.show()

plt.plot(error_list)
plt.show()

CPL.finalize()
MPI.Finalize()




