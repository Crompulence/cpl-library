import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

ppdir = '/home/asufian/Programs/coupled_LAMMPS_OpenFOAM/pyDataView'
sys.path.append(ppdir)
import postproclib as ppl

#Get Post Proc Object
fdir = './openfoam/'
PPObj = ppl.OpenFOAM_PostProc(fdir)

#Get plotting object
pObj = PPObj.plotlist['p']
UbObj = PPObj.plotlist['Ub']

#Get profile
h, p = pObj.profile(axis=1, startrec=pObj.maxrec, endrec=pObj.maxrec)
h, Ub = UbObj.profile(axis=1, startrec=UbObj.maxrec, endrec=UbObj.maxrec)
Ub = Ub[:,1]

#Sim box sizes
nx = pObj.Raw.ncx
ny = pObj.Raw.ncy
nz = pObj.Raw.ncz
xL = pObj.Raw.xL
yL = pObj.Raw.yL
zL = pObj.Raw.zL

#Simulation parameters
data = np.loadtxt('lammps/fcc0.dump', skiprows=9)
dp = np.mean(data[:, 10]);

porousStart = 91
porousEnd = 108

phi = 0.74
rho = 2.65
mu = 1.e-2
K = 435.
Ubb = 0.01

rp = 0.5*dp
eps = 1 - phi
L = (porousEnd - porousStart - 1)*(yL/ny)
cvol = (xL/nx)*(yL/ny)*(zL/nz)
pvol = cvol*phi
cnp = pvol/((np.pi/6)*dp**3)
Cd = 3*np.pi*mu*dp*K
cCd = cnp*Cd
pBot = 4.5*mu*phi*K*L*Ubb/(rp**2)

#Velocity plot
UbSol = np.concatenate(
    (Ubb*np.ones([porousStart, 1]), 
    (Ubb/eps)*np.ones([porousEnd-porousStart, 1]), 
    Ubb*np.ones([ny-porousEnd, 1])))

plt.plot(h, Ub, 'r-')
plt.plot(h, UbSol, 'k--')
plt.xlabel('Height (cm)')
plt.ylabel('Velocity (cm/s)')
plt.legend(('Numerical', 'Expected Solution'))
plt.savefig('fig_velocityfield.png')
plt.close()

#Pressure plot
pSol = np.concatenate(
    (pBot*np.ones([porousStart]), 
    np.linspace(pBot, 0, porousEnd-porousStart), 
    np.zeros([ny-porousEnd])))

#Plot only normal component
plt.plot(h, p, 'r-')
plt.plot(h, pSol, 'k--')
plt.xlabel('Height (cm)')
plt.ylabel('Pressure (Ba = 0.1Pa)')
plt.legend(('Numerical', 'Expected Solution'))
plt.savefig('fig_pressurefield.png')
