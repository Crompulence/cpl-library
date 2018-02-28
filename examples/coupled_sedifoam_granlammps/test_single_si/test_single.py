import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import os

OPEN_FOAM_CASE = 'openfoam/'
LAMMPS_CASE = 'lammps/'
OPEN_FOAM_BIN = '../CPLSediFOAM'
LAMMPS_BIN = '../lmp_cpl'
LAMMPS_INPUT = 'single.in'

os.system('./run.sh -L ' + LAMMPS_BIN + ' -O ' + OPEN_FOAM_BIN + ' -i ' + LAMMPS_INPUT)

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

#Load thermo_output data into numpy array
fName = LAMMPS_CASE + 'thermo_output.txt'
data = np.loadtxt(fName, skiprows=1)
t = data[:, 0]*dt
xp = data[:, 1:4]
xp = xp[:, gIdx]
vp = data[:, 4:7]
vp = vp[:, gIdx]

#Analytical Stokes solution
mp = rhop*(math.pi/6)*dp**3
mf = rhof*(math.pi/6)*dp**3
kd = 3*math.pi*mu*dp

vel = ((mp - mf)*g/kd)*(1 - np.exp(-kd*t/mp)) + v0*np.exp(-kd*t/mp)
disp = ((mp - mf)*g/kd)*(t - (mp/kd)*(1 - np.exp(-kd*t/mp))) + (v0*mp/kd)*(1 - np.exp(-kd*t/mp)) + x0

#Compare numerical and analytical solution
tol = 0.01
errvel = abs((vel - vp)/vel) <= tol
errdisp = abs((disp - xp)/disp) <= tol

def test_displacement():
	idx = np.where(errdisp == False)[0]
	print('Displacement profile at time t = ' + str(t[idx]) + ' exceeds the specified error tolerance of ' + str(tol*100) + '%')
	assert(errdisp.all())

def test_velocity():
    idx = np.where(errvel == False)[0]
    print('Velocity profile at time t = ' + str(t[idx]) + ' exceeds the specified error tolerance of ' + str(tol*100) + '%')
    assert(errvel.all())

def test_terminal():
    print('Terminal velocity not attained (within ' + str(tol*100) + '% tolerance). Velocity at end of simulation is ' + str(vp[-1]) + ', while terminal velocity is ' + str(((mp - mf)*g/kd)) + '.')
    assert(abs((vp[-1] - ((mp - mf)*g/kd))/((mp - mf)*g/kd) <= tol))

#Plot comparison of velocity and displacement profile
plt.plot(np.insert(t, 0, 0), np.insert(vp, 0, 0), 'r-', linewidth=3.0)
plt.plot(np.insert(t, 0, 0), np.insert(vel, 0, 0), 'k-')
plt.plot(np.insert(t, 0, 0), ((mp - mf)*g/kd)*np.ones_like(np.insert(t, 0, 0)), 'k--')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.legend(('Numerical', 'Analytical', 'Terminal'))
plt.savefig('fig_velocity.png')
plt.close()

plt.plot(np.insert(t, 0, 0), np.insert(xp, 0, x0), 'r-', linewidth=3.0)
plt.plot(np.insert(t, 0, 0), np.insert(disp, 0, x0), 'k-')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend(('Numerical', 'Analytical'))
plt.savefig('fig_displacement.png')
plt.close()