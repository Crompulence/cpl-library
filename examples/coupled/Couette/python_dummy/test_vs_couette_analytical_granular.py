#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from CouetteAnalytical import CouetteAnalytical as CA

#Import CPL library
from cplpy import CPL

import sys
sys.path.append("/home/ubuntu/codes/pyDataView")
import postproclib as ppl

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

#Plot output
fig, ax = plt.subplots(1,1)
plt.ion()
plt.show()

#Analytical solution
dt = 0.05
U = 1.
nu = 1e-02
Re = xyzL[1]*U/nu
ncy = CPL.get("ncy")
dy = CPL.get("dy")
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)
y = np.linspace(dy/2., xyzL[1] - dy/2., num=ncy)

#Porous region solution
pcell = 4
yl_cfd = CPL.get("yl_cfd")
Lmin = (pcell+0.5)*dy
Re = (xyzL[1]-Lmin)*U/nu
CApObj = CA(Re=Re, U=U, Lmin=Lmin, Lmax=xyzL[1], npoints=ncy)

# Coupling parameters
eps = 0.0001
A = 1./eps
l2 = []

#Main time loop
pcount = 0
for time in range(1000):

    print(time)

    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e]
    sendbuf[7,:,:,:] = 0.0

    #Set pcell to zero porosity -- this needs to be two cells
    #sendbuf[7,:,pcell:pcell+2,:] = 1.0 - eps

    # Set porous force based on drag difference between velocity
    sendbuf[3,:,pcell:pcell+1,:] = 0.05*( sendbuf[0,:,pcell:pcell+1,:]
                                         -recvbuf[0,:,pcell:pcell+1,:])
    CPL.send(sendbuf)

    #Plot data
    if time%20 == 0:

        rec = int(time/20)

        try:
            #Plot values from OpenFOAM read by postproc library 
            OpenFOAMuObj = ppl.OpenFOAM_vField("./openfoam", parallel_run=True)
            y, u = OpenFOAMuObj.profile(1, startrec=rec, endrec=rec)
            l, = ax.plot(u[:,0], y, 'ro-', label="OpenFOAM domain from file")

            #Plot values recieved from SediFOAM
            ur = np.mean(recvbuf[0,:,:,:],(0,2))
            ax.plot(ur, y, 'g^', ms=10, label="Recieved value from SediFOAM")

            #Plot analytical solution
            y_anal, u_anal = CAObj.get_vprofile(time*dt)
            ax.plot(u_anal[:], y_anal, 'k.-', label="Analytical Solution")

            y_p, u_p = CApObj.get_vprofile(time*dt)
            ax.plot(u_p, y_p, 'kx--', label="Analytical Solution porous")

            ax.set_xlim([-0.1,1.1])
            ax.set_xlabel("$u$",fontsize=24)
            ax.set_ylabel("$y$",fontsize=24)

            plt.legend(loc=4)
            plt.pause(0.1)

            #plt.savefig('out{:05}.png'.format(pcount));
            pcount += 1

            ax.cla()
            #print("MDTime=", time, OpenFOAMuObj.maxrec, rec)

        except ppl.pplexceptions.DataNotAvailable:
            print("Data not found", time, rec)

        except ppl.pplexceptions.OutsideRecRange:
            print("Outside record range", time, rec, OpenFOAMuObj.maxrec)

CPL.finalize()
MPI.Finalize()




