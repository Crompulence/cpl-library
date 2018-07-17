import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import os

fDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(fDir)
OPEN_FOAM_CASE = 'openfoam/'
LAMMPS_CASE = 'lammps/'
os.system('./run.sh')

#Extract parameters used in the simulations: timestep (dt), gravity (g), direction of gravity (gIdx),
#density of fluid (rhof), dynamic viscosity of fluid (mu), particle diameter (dp), density of particle (rhop),
#initial position (x0) and velocity (v0) in the direction of gravity
with open(OPEN_FOAM_CASE + 'constant/environmentalProperties', 'r') as fObj:
    for l in fObj:
        if l.startswith('g'):
            g = l[l.find('(')+1:l.find(')')]
            g = [float(i) for i in g.split()]
            gIdx = [i for i, x in enumerate(g) if (x != 0)][0]
            g = g[gIdx]

with open(OPEN_FOAM_CASE + 'constant/transportProperties', 'r') as fObj:
    for l in fObj:
        if 'rhob' in l:
            rhof = l.split()[-1]
            rhof = float(rhof[:-1])
        if 'nub' in l:
            nu = l.split()[-1]
            nu = float(nu[:-1])
    mu = nu*rhof

with open(LAMMPS_CASE + 'single.lj', 'r') as fObj:
    lines = filter(None, (line.rstrip() for line in fObj))   
    for i, l in enumerate(lines):
        if l == 'Atoms':
            data = lines[i+1].split()
            dp = float(data[2])
            rhop = float(data[3])
            x0 = data[4:7]
            x0 = float(x0[gIdx])
            
        if l == 'Velocities':
            data = lines[i+1].split()
            v0 = data[1:4]
            v0 = float(v0[gIdx])

with open(LAMMPS_CASE + 'single.in', 'r') as fObj:
    for l in fObj:
        if 'timestep' in l:
        	dt = float(l.split()[1])
        if 'pair_style' in l:
            data = l.split()
            Gp = float(data[-3])
            nup = float(data[-2])
        if 'wall/gran' in l:
            data = l.split()
            Gw = float(data[5])
            nuw = float(data[6])

#Load thermo_output data into numpy array
fName = LAMMPS_CASE + 'thermo_output.txt'
data = np.loadtxt(fName, skiprows=1)
t = data[:, 0]*dt
xp = data[:, 1:4]
xp = xp[:, gIdx]
vp = data[:, 4:7]
vp = vp[:, gIdx]

#Analytical solution for equilibrium position
mp = rhop*(math.pi/6)*dp**3
mf = rhof*(math.pi/6)*dp**3
kb = (mp - mf)*g

rp = 0.5*dp
Ep = 2*Gp*(1 + nup)
Ew = 2*Gw*(1 + nuw)
Estar = 1/(((1 - nup**2)/Ep) + ((1 - nuw**2)/Ew))
kc = (4./3.)*Estar*np.sqrt(rp)

delta = np.power((-kb/kc), 2./3.)
yequil = rp - delta

#Compare numerical and analytical solution
tol = 0.01
def test_equilibrium():
    print('Equilibrium position (' + str(xp[-1]) + ') did not match analytical solution (' + str(yequil) + ') within tolerance of ' + str(100*tol) + '%')
    assert(abs((xp[-1] - yequil)/yequil) < tol)

#Plot comparison of velocity and displacement profile
plt.plot(np.insert(t, 0, 0), np.insert(vp, 0, 0), 'r-', linewidth=3.0)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (cm/s)')
plt.legend(('Numerical', 'Terminal Velocity'))
plt.savefig('fig_velocity.png')
plt.close()

plt.plot(np.insert(t, 0, 0), np.insert(xp, 0, x0), 'r-', linewidth=3.0)
plt.plot(np.insert(t, 0, 0), yequil*np.ones_like(np.insert(t, 0, 0)), 'k--')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (cm)')
plt.legend(('Numerical', 'Analytical'))
plt.savefig('fig_displacement.png')
plt.close()