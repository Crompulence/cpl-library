
import matplotlib.pyplot as plt
import numpy as np
import sys
from CouetteAnalytical import CouetteAnalytical as CA

ppdir = '/home/es205/codes/python/pyDataView/'
sys.path.append(ppdir)
import postproclib as ppl

# Font size and type
fsize = 24
plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#Get Post Proc Object
fdir = './'
PPObj = ppl.CPL_PostProc(fdir)
print(PPObj)

#Get plotting object
plotObj = PPObj.plotlist['CPL Velocity']
CFDObj = plotObj.cfd_field
MDObj = plotObj.md_field

#Get domain sizes and setup couette analytical soln
dx, dy, dz = plotObj.cpl_dxyz
X, Y, Z = plotObj.get_cpl_grid()
xyzL = [X.max(), Y.max(), Z.max()]
nxyz = [X.shape[0], Y.shape[0], Z.shape[0]]
U = 1.0
nu = CFDObj.Raw.nu[0]
Re = (xyzL[1]-dy)/nu
rec_dt = float(CFDObj.Raw.header.headerDict['controlDict']['writeInterval'])
 
fig, ax = plt.subplots(1,1)
#plt.ion()
#plt.show()

ave = 0
for rec in [25,50,100,250,790]:#plotObj.maxrec-1-ave]:

    #Get coupled solutions
    yMD, uMD, yCFD, uCFD = plotObj.profile_both(axis=1, startrec=rec-ave, endrec=rec+ave)

    ymin = yMD[0]
#    ax.plot(uMD[:-12,0], yMD[:-12]-4.33021, 'ro')
#    ax.plot(uCFD[:,0], yCFD-9.296, 'bx')

#    CAObj = CA(Re=Re, U=U, Lmin=0.0, Lmax=yCFD[-1]-9.296+5.66/2., npoints=2*nxyz[1]+2) 
#    y_anal, u_anal = CAObj.get_vprofile(rec*rec_dt)
#    ax.plot(u_anal, y_anal, 'k-')

    #yMD = np.linspace(-4, 50, uMD.shape[0])
    ax.plot(uMD[:,0], yMD[:], 'bo')
    ax.plot(uCFD[:,0], yCFD, 'gx')

    haloCFD = CFDObj.read_halo(startrec=rec,endrec=rec, 
                               haloname="movingWall")
    haloCPL = CFDObj.read_halo(startrec=rec, endrec=rec, 
                               haloname="CPLReceiveMD")
    ax.plot(np.mean(haloCPL[...,0],(0,2,3)), yCFD[0 ]-dy/2., 'gs')
    ax.plot(np.mean(haloCFD[...,0],(0,2,3)), yCFD[-1]+dy/2., 'gs')

    ax.plot(uMD[6,0], yMD[6], 'rs')

    #Get analytical solns
    CAObj = CA(Re=Re, U=U, Lmin=ymin+2., Lmax=yCFD[-1]+dy/2., npoints=2*nxyz[1]+2)
    y_anal, u_anal = CAObj.get_vprofile(rec*rec_dt)
    ax.plot(u_anal, y_anal, 'k-')

    
    #plt.pause(0.001)
    #plt.cla()
plt.ylabel("$y$")
plt.xlabel("$u$")
plt.show()
#plt.savefig("Coupled_Couette.eps", bbox_inches="tight")
#plt.savefig("Coupled_Couette.pdf", bbox_inches="tight")
