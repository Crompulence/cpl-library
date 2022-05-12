import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

#Particle
D = 3.5e-4
r = D/2.
V=(4./3.)*np.pi*r**3.
rho = 2.65e3
E1 = 2.9e4
nu1=0.2
f1 = 0.25
mass = rho*V

#Wall
E2 = E1
nu2 = nu1

#Interaction
E = 1./( (1.-nu1**2)/E1 + (1.-nu2**2)/E2)
k = (4./3.)*E*r**0.5


skip = 20
def heaviside(x):
    return 1 * (x > 0)

def analytical_gravity(z0, v0, t0, t):
    g = 9.81
    return z0 + v0*(t-t0) - 0.5*g*(t-t0)**2

def analytical_interaction(k, m, v0, t):
    ta = t - t[0]
    omega = np.sqrt(k/m)
    return 0.5*v0*np.sin(omega*ta)

def get_maxdisp_during_bounce(vin):

    return ((3/8.)*(mass*vin**2)/(E*r**0.5))**(2/5.)


#Get data from file
with open('./log.lammps') as f:
    for l in f.readlines():
        if l.find("timestep") != -1:
            dt = float(l.strip('timestep'))
            break
    
print(("timestep = ", dt))
data = np.genfromtxt('./thermo_output.txt')
t = data[:,0]*dt
z = data[:,1]
v = data[:,2]
f = data[:,3]
plt.plot(t[::-skip], z[::-skip], 'ro',ms=10.)

#Analytical solutions
z0 = z.max()

#tm = (2.*(z0-D)/g)**.5
#vr = g*tm
zp = np.copy(z)
zp[z<D/2.]=1000.
mins = argrelextrema(zp, np.less)[0]

#Plot minimum points
for m in mins:
    plt.plot(t[m], z[m], 'ko', ms=10)
    maxdisp = get_maxdisp_during_bounce(v[m])
    plt.plot(t[m], z[m]-maxdisp, 'ko', ms=10)
    print(maxdisp)

m = mins[0]
mp1 = mins[1]
ta = t[m:mp1]
analytical_interaction(k, mass, v[m], ta)



#Plot gravity bit
error = []
for i in range(1,len(mins)-1,2):
    m = mins[i]
    mp1 = mins[i+1]
    ta = t[m:mp1]
    za = analytical_gravity(z[m], v[m], t[m], ta)
    plt.plot(ta, za, 'b-', lw=3.)
    error.append([ta, z[m:mp1]-za])

plt.show()

