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
        self.v = np.random.randn(self.v.shape[0], self.v.shape[1])
        vsum = np.sum(self.v,0)/self.N
        for i in range(self.N):
            self.v[i,:] -= vsum

        #Setup velocity averaging
        self.xbin = 8; self.ybin = 8
        self.xb = np.linspace(-self.halfdomain[0],
                               self.halfdomain[0], self.xbin)
        self.yb = np.linspace(-self.halfdomain[1],
                               self.halfdomain[1], self.ybin)
        self.Xb, self.Yb = np.meshgrid(self.xb, self.yb)

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

    def get_velfield(self, bins):

        velbin = np.zeros([bins[0],bins[1],2])
        binsize = self.domain/bins
        for i in range(self.r.shape[0]):
            ib = [int((self.r[i,ixyz]+0.5*self.domain[ixyz])
                       /binsize[ixyz]) for ixyz in range(self.nd)]
            velbin[ib[0],ib[1],:] += self.v[i,:]

        return velbin

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
                ib = [int((self.r[i,ixyz]+self.halfdomain[ixyz])
                           /cellsize[ixyz]) for ixyz in range(nd)]
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

            #Peridic boundary
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


    def constraint_force(self, u_CFD, constraint_cell, alpha=0.1):

        #Get the MD velocity field
        binsize = self.domain/u_CFD.shape[0:2]
        binsize_MD = self.domain/[self.xbin,self.ybin]

        u_MD = self.get_velfield([self.xbin,self.ybin])
        
        #Extract CFD value
        F = np.zeros(2)
        for i in range(self.N):
            ib = [int((self.r[i,ixyz]+self.halfdomain[ixyz])
                       /binsize[ixyz]) for ixyz in range(self.nd)]
            if ib[0] > u_CFD.shape[1]:
                ib[0] = u_CFD.shape[1]
            if ib[1] > u_CFD.shape[2]:
                ib[1] = u_CFD.shape[2]
            F[:] = alpha*(u_CFD[ib[0],ib[1],:] - u_MD[ib[0],ib[1],:])
            self.a[i,:] += F[:]

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
        draw_grid(ax, nx=self.xbin+1, ny=self.ybin+1, nz=1,
                  xmin=-self.halfdomain[0], xmax=self.halfdomain[0],
                  ymin=-self.halfdomain[1], ymax=self.halfdomain[1])

        #Get velocity field
        self.velbin = self.get_velfield([self.xbin,self.ybin])

        #Plot velocity profile offset to the left
        axisloc = self.halfdomain[0]+1
        ax.arrow(axisloc,-self.halfdomain[1], 0.,self.domain[1],  
                 width=0.015, color="k", clip_on=False, head_width=0.12, head_length=0.12)
        ax.arrow(axisloc-1,0., 2.,0.,  width=0.015, 
                 color="k", clip_on=False, head_width=0.12, head_length=0.12)
        ax.plot(np.mean(self.velbin[:,:,0],0)+axisloc,self.yb,'g-')

    #    cb=ax.imshow(velbin[:,:,0],interpolation="none",
    #                 extent=[-halfdomain[0],halfdomain[0],
    #                         -halfdomain[1],halfdomain[1]], 
    #                cmap=plt.cm.RdYlBu_r,vmin=-3.,vmax=3.)
    #    if first_time:
    #        plt.colorbar(cb)
    #        first_time=False

        if showarrows:
            #Show velocity of molecules
            Q = ax.quiver(self.r[:,0], self.r[:,1],
                          self.v[:,0], self.v[:,1], color='k')

        #Set limits and plot
        ax.set_xlim((-self.halfdomain[0], self.halfdomain[0]+2.))
        ax.set_ylim((-self.halfdomain[1], self.halfdomain[1]))
        plt.pause(0.001)
        plt.cla()


if __name__ == "__main__":

    Nsteps = 10000

    md = MD()

    #Main run
    for step in range(Nsteps):

        print("MD time = ", step, " of ", Nsteps)

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



