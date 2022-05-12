import numpy as np
import matplotlib.pyplot as plt

from draw_grid import draw_grid



class MD:

    def __init__(self,
                 initialnunits = [3,8],
                 density = 0.8,
                 nd = 2,
                 rcutoff = 2.**(1./6.),
                 dt = 0.005,
                 Tset = 1.3,               #After equilbirum approx temp of 1
                 forcecalc = "allpairs",
                 wallwidth = [0.,0.],
                 wallslide = [0.,0.],
                 newfig=None):

        self.initialnunits = initialnunits
        self.density = density
        self.nd = nd
        self.rcutoff = rcutoff
        self.dt = dt
        self.forcecalc = forcecalc
        self.wallwidth = np.array(wallwidth)
        self.wallslide = np.array(wallslide)

        self.rcutoff2 = rcutoff**2
        self.first_time=True
        self.tstep = 0
        self.time = 0.
        self.periodic = [True, True]
        self.spec_wall = [False, False]

        #Setup Initial crystal
        self.domain = np.zeros(2)
        self.volume=1.    #Set domain size to unity for loop below
        for ixyz in range(nd):
            self.domain[ixyz] = (initialnunits[ixyz]/((density/4.0)**(1.0/nd)))
            self.volume = self.volume*self.domain[ixyz]        #Volume based on size of domain
        self.halfdomain = 0.5*self.domain

        #Allocate arrays
        self.N = 4*initialnunits[0]*initialnunits[1]
        self.tag = np.zeros(self.N, dtype=int)
        self.r = np.zeros((self.N,2))
        self.v = np.zeros((self.N,2))
        self.a = np.zeros((self.N,2))
        self.v = Tset*np.random.randn(self.v.shape[0], self.v.shape[1])
        vsum = np.sum(self.v,0)/self.N
        for i in range(self.N):
            self.v[i,:] -= vsum

        #Setup velocity averaging
        self.veluptodate = 0
        self.xbin = 8; self.ybin = 8
        self.dx = self.domain[0]/self.xbin
        self.dy = self.domain[1]/self.ybin
        self.xb = np.linspace(-self.halfdomain[0],
                               self.halfdomain[0], self.xbin)
        self.yb = np.linspace(-self.halfdomain[1],
                               self.halfdomain[1], self.ybin)
        self.Xb, self.Yb = np.meshgrid(self.xb, self.yb)
        self.mbin = np.zeros([self.xbin, self.ybin])
        self.velbin = np.zeros([2, self.xbin, self.ybin])

        self.setup_crystal()
        self.setup_walls(wallwidth)

        if newfig == None:
            self.fig, self.ax = plt.subplots(1,1)
            plt.ion()
            plt.show()

    #Molecules per unit FCC structure (3D)
    def setup_crystal(self):
        initialunitsize = self.domain / self.initialnunits
        n  = 0      #Initialise global N counter n
        c = np.zeros(2); rc = np.zeros(2)
        for nx in range(1,self.initialnunits[0]+1):
            c[0] = (nx - 0.750)*initialunitsize[0]
            for ny in range(1,self.initialnunits[1]+1):
                c[1] = (ny - 0.750)*initialunitsize[1] 
                for j in range(4):    #4 Molecules per cell

                    rc[:] = c[:]
                    if j is 1:
                        rc[0] = c[0] + 0.5*initialunitsize[0]
                    elif j is 2:        
                        rc[1] = c[1] + 0.5*initialunitsize[1]
                    elif j is 3:
                        rc[0] = c[0] + 0.5*initialunitsize[0]
                        rc[1] = c[1] + 0.5*initialunitsize[1]
                
                    #Correct to local coordinates
                    self.r[n,0] = rc[0]-self.halfdomain[0]
                    self.r[n,1] = rc[1]-self.halfdomain[1]
                    n += 1   #Move to next particle

    def setup_walls(self, wallwidth):

        for i in range(self.N):
            if (self.r[i,1]+self.halfdomain[1] < wallwidth[0]):
                self.tag[i] = 1
            elif (self.r[i,1]+self.halfdomain[1] > self.domain[1]-wallwidth[1]):
                self.tag[i] = 1
            else:
                self.tag[i] = 0

        if any(self.wallwidth > 0.):
            self.periodic[1] = False
            self.spec_wall[1] = True

    def LJ_accij(self,rij2):

        invrij2 = 1./rij2
        return 48.*(invrij2**7-.5*invrij2**4)


    def get_bin(self, r, binsize):
        ib = [int((r[ixyz]+0.5*self.domain[ixyz])
              /binsize[ixyz]) for ixyz in range(self.nd)]
        return ib

    def get_velfield(self, bins, freq=25, plusdt=False, getmbin=False):

        #Update velocity if timestep dictates 
        if ((self.tstep%freq == 0)
             and (self.tstep != self.veluptodate)):
            mbin = np.zeros([bins[0], bins[1]])
            velbin = np.zeros([2, bins[0], bins[1]])
            binsize = self.domain/bins
            #Loop over all molecules in r
            for i in range(self.r.shape[0]):
                ib = self.get_bin(self.r[i,:], binsize)
                mbin[ib[0], ib[1]] += 1
                if plusdt:
                    vi = self.v[i,:] + self.dt*self.a[i,:]
                else:
                    vi = self.v[i,:]
                velbin[:, ib[0], ib[1]] += vi
            self.mbin = mbin
            self.velbin = velbin
            self.veluptodate = self.tstep
        else:
            mbin = self.mbin
            velbin = self.velbin

        u = np.divide(velbin,mbin) 
        u[np.isnan(u)] = 0.0
        if getmbin:
            return u, mbin
        else:
            return u

    def force(self, showarrows=False, ax=None):

        if ax == None:
            ax=self.ax

        #Force calculation
        self.a = np.zeros((self.N,2))
        if self.forcecalc is "allpairs":

            for i in range(self.N):
                for j in range(i+1,self.N):

                    ri = self.r[i,:]; rj = self.r[j,:]
                    rij = ri - rj

                    #Nearest neighbour
                    for ixyz in range(self.nd):
                        if self.periodic[ixyz]:
                            if (np.abs(rij[ixyz]) > self.halfdomain[ixyz]):
                                rij[ixyz] -= np.copysign(self.domain[ixyz],rij[ixyz]) 

                    #Get forces
                    rij2 = np.dot(rij,rij)
                    if rij2 < self.rcutoff2:
                        fij = self.LJ_accij(rij2)*rij
                        self.a[i,:] += fij
                        self.a[j,:] -= fij

                        if showarrows:
                            ax.quiver(self.r[i,0],self.r[i,1],
                                      rij[0], rij[1],color='red',width=0.002)
                            ax.quiver(self.r[j,0],self.r[j,1], 
                                     -rij[0],-rij[1],color='red',width=0.002)

        elif "celllist":

            #Build celllist
            self.cells = [int(self.domain[i]/(self.rcutoff)) for i in range(self.nd)]
            self.cell = np.zeros(self.cells, dtype=object)
            for icell in range(self.cells[0]):
                for jcell in range(self.cells[1]):
                    self.cell[icell,jcell] = []

            cellsize = self.domain/cells
            for i in range(self.N):
                ib = self.get_bin(r[i,:], cellsize)
                self.cell[ib[0],ib[1]].append(i)

            #For each cell, check molecule i
            for icell in range(self.cells[0]):
                for jcell in range(self.cells[1]):
                    for i in cell[icell,jcell]:
                        ri = self.r[i,:]

                        #Check i against all molecules in adjacent cells
                        for aicell in [icell-1,icell,(icell+1)%cells[0]]:
                            for ajcell in [jcell-1,jcell,(jcell+1)%cells[1]]:
                                for j in cell[aicell,ajcell]:
                                    rj = self.r[j,:]
                                    rij = ri - rj

                                    #Nearest neighbour
                                    for ixyz in range(nd):
                                        if self.periodic[ixyz]:
                                            if (np.abs(rij[ixyz]) > self.halfdomain[ixyz]):
                                                rij[ixyz] -= np.copysign(self.domain[ixyz],rij[ixyz]) 

                                    #Get forces
                                    rij2 = np.dot(rij,rij)
                                    if rij2 < self.rcutoff2 and rij2 > 1e-8:
                                        fij = self.LJ_accij(rij2)*rij
                                        self.a[i,:] += fij

                                        if showarrows:
                                            ax.quiver(self.r[i,0],self.r[i,1],
                                                      -rij[0],-rij[1],color='red',width=0.002)


    def verlet(self):

        #Verlet time advance
        for i in range(self.N):

            if self.tag[i] == 0:
                self.v[i,:] += self.dt*self.a[i,:]
                self.r[i,:] += self.dt*self.v[i,:]
            else:
                #Fixed molecule, add sliding
                self.v[i,:] = self.wallslide[:]
                self.r[i,:] += self.dt*self.wallslide[:]

            #Peridic boundary and specular wall
            for ixyz in range(self.nd):
                if self.r[i,ixyz] > self.halfdomain[ixyz]:
                    if self.spec_wall[ixyz]:
                        overshoot = self.r[i,ixyz]-self.halfdomain[ixyz]
                        self.r[i,ixyz] -= 2.*overshoot
                        self.v[i,ixyz] = -self.v[i,ixyz] 

                    else:
                        self.r[i,ixyz] -= self.domain[ixyz]

                elif self.r[i,ixyz] < -self.halfdomain[ixyz]:
                    if self.spec_wall[ixyz]:
                        overshoot = -self.halfdomain[ixyz]-self.r[i,ixyz]
                        self.r[i,ixyz] += 2.*overshoot
                        self.v[i,ixyz] = -self.v[i,ixyz] 
                    else:
                        self.r[i,ixyz] += self.domain[ixyz]

        #Increment current time step
        self.tstep += 1
        self.time = self.tstep*self.dt 

    def constraint_force(self, u_CFD, constraint_cell, alpha=0.1):

        #Get the MD velocity field
        binsize_CFD = self.domain/u_CFD.shape[1:2]
        binsize_MD = self.domain/[self.xbin, self.ybin]
        assert binsize_CFD[0] == binsize_MD[0]
        assert binsize_CFD[1] == binsize_MD[1]
        u_MD, mbin = self.get_velfield([self.xbin,self.ybin], getmbin=True)
       
        #Extract CFD value
        F = np.zeros(2)
        ucheck = np.zeros([2,self.xbin])
        hd = self.halfdomain
        for i in range(self.N):
            ib = self.get_bin(self.r[i,:], binsize_MD)
            #Ensure within domain
            if ib[0] > u_MD.shape[1]:
                ib[0] = u_MD.shape[1]
            if ib[1] > u_MD.shape[2]:
                ib[1] = u_MD.shape[2]
            #only apply to constrained cell
            if ib[1] == constraint_cell:
                F[:] = alpha*(u_CFD[:,ib[0],0] - u_MD[:,ib[0],ib[1]])
                if (mbin[ib[0],ib[1]] != 0):
                    self.a[i,:] += F[:]/float(mbin[ib[0],ib[1]])
                else:
                    pass
                self.ax.quiver((ib[0]+.5)*self.dx-hd[0],
                             (ib[1]+.5)*self.dy-hd[1],F[0],F[1],
                              color='red',angles='xy',scale_units='xy',scale=1)



    def CV_constraint_force(self, u_CFD, constraint_cell):

        #Get the MD velocity field
        binsize_CFD = self.domain/u_CFD.shape[1:2]
        binsize_MD = self.domain/[self.xbin, self.ybin]
        assert binsize_CFD[0] == binsize_MD[0]
        assert binsize_CFD[1] == binsize_MD[1]
        u_MD = self.get_velfield([self.xbin,self.ybin], freq=1)
       
        #Apply force value
        F = np.zeros(2)
        du_MDdt = np.zeros(2)
        du_CFDdt = np.zeros(2)
        hd = self.halfdomain
        for i in range(self.N):
            ib = self.get_bin(r[i,:], binsize_MD)
            #Ensure within domain
            if ib[0] > u_MD.shape[1]:
                ib[0] = u_MD.shape[1]
            if ib[1] > u_MD.shape[2]:
                ib[1] = u_MD.shape[2]
            #only apply to constrained cell
            if ib[1] == constraint_cell:
                du_MDdt[:] = (u_MD[:,ib[0],ib[1]] 
                       - self.u_MD[:,ib[0],ib[1]])
                du_CFDdt[:] = (u_CFD[:,ib[0],ib[1]] 
                        - self.u_CFD[:,ib[0],ib[1]])

                if (mbin[ib[0],ib[1]] != 0):
                    F[:] = ( (du_MDdt - du_CFDdt)
                           /(float(mbin[ib[0],ib[1]])*self.dt))
                else:
                    F[:] = 0.

        self.u_MD = u_MD
        self.u_CFD = u_CFD

    #Plot molecules
    def plot(self, ax=None, showarrows=False):

        if ax == None:
            ax=self.ax

        for i in range(self.N):
            if (self.tag[i] == 0):
                ax.plot(self.r[i,0],self.r[i,1],'ko',alpha=0.5)
            else:
                ax.plot(self.r[i,0],self.r[i,1],'ro', ms=7.)

        #Overlay grid
        draw_grid(ax, nx=self.xbin, ny=self.ybin, nz=1,
                  xmin=-self.halfdomain[0], xmax=self.halfdomain[0],
                  ymin=-self.halfdomain[1], ymax=self.halfdomain[1])

        #Get velocity field
        u = self.get_velfield([self.xbin,self.ybin])

        #Plot velocity profile offset to the left
        axisloc = self.halfdomain[0]+1
        ax.arrow(axisloc,-self.halfdomain[1], 0.,self.domain[1],  
                 width=0.015, color="k", clip_on=False, head_width=0.12, head_length=0.12)
        ax.arrow(axisloc-1,0., 2.,0.,  width=0.015, 
                 color="k", clip_on=False, head_width=0.12, head_length=0.12)

        yp = np.linspace(-self.halfdomain[1]+.5*self.dy, self.halfdomain[1] - 0.5*self.dy, self.ybin)
        ax.plot(np.mean(u[0,:,:],1)+axisloc,yp,'g-x')

        sm = ax.imshow(u[0,:,:].T,aspect='auto',origin='lower',
                       extent=[-self.halfdomain[0], self.halfdomain[0],
                               -self.halfdomain[1], self.halfdomain[1]],
                       interpolation="none",vmin=-1.,vmax=1.,
                       alpha=0.5, cmap=plt.cm.RdYlBu_r)

#        sm = ax.pcolormesh(self.Xb,self.Yb,u[0,:,:].T,vmin=-1.,vmax=1.,alpha=0.5,
#                          cmap=plt.cm.RdYlBu_r)

#        cb=ax.imshow(u[0,:,:],interpolation="none",
#                     extent=[-self.halfdomain[0],self.halfdomain[0],
#                             -self.halfdomain[1],self.halfdomain[1]], 
#                    cmap=plt.cm.RdYlBu_r,vmin=-3.,vmax=3.)
        if self.first_time:
            plt.colorbar(sm)
            self.first_time=False

        if showarrows:
            #Show velocity of molecules
            Q = ax.quiver(self.r[:,0], self.r[:,1],
                          self.v[:,0], self.v[:,1], color='k')

        #Set limits and plot
        ax.set_xlim((-self.halfdomain[0], self.halfdomain[0]+2.))
        ax.set_ylim((-self.halfdomain[1], self.halfdomain[1]))
        plt.pause(0.001)
        plt.cla()

        print(("Temperature =", np.sum(self.v[:,0]**2+self.v[:,1]**2)/(2.*self.N)))


if __name__ == "__main__":

    Nsteps = 10000

    md = MD()

    #Main run
    for step in range(Nsteps):

        print(("MD time = ", md.tstep, " of ", Nsteps))

        md.force()

        #=======================================================
        # Call to CPL-LIBRARY goes here to
        # recieve u_CFD in constraint region
        # and force is applied
        # F = (1/tau)*(u_CFD - u_MD)
        #=======================================================

        md.verlet()

        #=======================================================
        #Call to CPL-LIBRARY goes here to send u_MD at boundary
        #=======================================================

        md.plot()



