import numpy as np
from mpi4py import MPI
import sys

#Import CPL library
from cplpy import CPL

#initialise MPI
comm = MPI.COMM_WORLD

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([2., 10., 2.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#initialise CPL
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)
CPL.setup_md(MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]]), xyzL, xyz_orig)

#Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=9, send_size=8)

#Main time loop
porousStart = 91
porousEnd = 108

phi = 0.74
dp = 0.05
rho = 2.65
mu = 1.e-2
K = 435.
Ubb = 0.01

rp = 0.5*dp
eps = 1 - phi
L = (porousEnd - porousStart - 1)*(xyzL[1]/np.shape(recvbuf)[2])
cvol = (xyzL[0]/np.shape(recvbuf)[1])*(xyzL[1]/np.shape(recvbuf)[2])*(xyzL[2]/np.shape(recvbuf)[3])
pvol = cvol*phi
cnp = pvol/((np.pi/6)*dp**3)
Cd = 3*np.pi*mu*dp*K
cCd = cnp*Cd
pBot = 4.5*mu*phi*K*L*Ubb/(rp**2)

print('Number of particles in CFD cell: ' + str(cnp))
print('Sum of drag coeff in CFD cell: ' + str(cCd))
print('Expected velocity through packed bed: ' + str(Ubb/eps))
print('Expected pressure gradient through packed bed: ' + str(pBot/L))
print('Expected pressure value on bottom boundary: ' + str(pBot))

for time in range(101):

    print(time)
    # Recv data: 
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    recvbuf, ierr = CPL.recv(recvbuf)

    # Zero send buffer and set porosity to one
    # [Ux, Uy, Uz, Fx, Fy, Fz, Cd, e] 
    sendbuf[:,:,:,:] = 0.

    sendbuf[6,:,porousStart:porousEnd,:] = cCd
    sendbuf[7,:,porousStart:porousEnd,:] = pvol

    CPL.send(sendbuf)

CPL.finalize()
MPI.Finalize()




