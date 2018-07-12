import matplotlib.pyplot as plt
import numpy as np
from CouetteAnalytical import CouetteAnalytical as CA

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([1., 100., 1.], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#Analytical solution
dt = 0.1
U = 1.
nu = 1.004e-2
Re = U/nu

ncy = 6
CAObj = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)
ncy = 20
CAObj2 = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)
ncy = 100
CAObj3 = CA(Re=Re, U=U, Lmin=0., Lmax=xyzL[1], npoints=2*ncy+1, nmodes=100*ncy)

#Main time loop
pcount = 0
for time in [1,10,100,500,1000,10000]:

    print(time)
    y_anal, u_anal = CAObj.get_vprofile(time*dt)
    plt.plot(u_anal,y_anal, 'k.-', label="Analytical Solution")
    y_anal, u_anal = CAObj2.get_vprofile(time*dt)
    plt.plot(u_anal,y_anal, 'ro-', label="Analytical Solution")
    y_anal, u_anal = CAObj3.get_vprofile(time*dt)
    plt.plot(u_anal,y_anal, 'bo-', label="Analytical Solution")

plt.show()

