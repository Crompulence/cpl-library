import numpy as np
import matplotlib.pyplot as plt

from draw_grid import draw_grid

class CFD:
    """
        Solve the diffusion equation
        du/dt = (rho/gamma) * d2u/dx2
    """
    def __init__(self, dt, nu = 1., 
                 xsize = 10, ysize = 10,
                 xmin = 0., xmax = 1.,
                 ymin = 0., ymax = 1.,
                 fig=None):

        #Define coefficients
        self.nu = nu
        self.xsize = xsize
        self.ysize = ysize
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.dt = dt

        #Define arrays
        self.x = np.linspace(xmin,xmax,xsize)
        self.dx = np.mean(np.diff(self.x))
        self.y = np.linspace(ymin,ymax,ysize)
        self.dy = np.mean(np.diff(self.y))

        #For plotting
        self.X,self.Y = np.meshgrid(self.x,self.y)

        #Check CFL stability conditions
        self.CFL =  (1./(2.*nu))*(self.dx*self.dy)**2/(self.dx**2+self.dy**2)
        if self.dt > self.CFL:
            print("Stability conditions violated, CFL=", self.CFL ,
                  "> dt=", self.dt," adjust dt, nu or grid spacing")
            quit()
        else:
            print("Timestep dt = ", self.dt, " CFL number= ", self.CFL)

        #initial condition
        self.u0 = np.zeros([xsize,ysize])

        #Setup first times
        self.u = np.copy(self.u0)
        self.u_mdt =  np.copy(self.u0)
        self.u_m2dt =  np.copy(self.u0)
        self.first_time = True

        #Setup figure
        if fig ==  None:
            self.fig, self.ax = plt.subplots(1,1)
            plt.ion()
            plt.show()
        else:
            self.fig = fig
            self.ax = fig.axes

    def set_bc(self, topwall=1., bottomwall=0.):
        #Enforce boundary conditions
        self.u[:,0] = bottomwall; self.u[:,-1] = topwall
        #Periodic boundaries
        self.u[-1,:] = self.u[1,:]; self.u[0,:] = self.u[-2,:]

    def update_time(self):
        # Save previous value
        self.u_m2dt = np.copy(self.u_mdt)
        self.u_mdt = np.copy(self.u)

        #Solve for new u
        for i in range(1,self.x.size-1):
            for j in range(1,self.y.size-1):
                #Diffusion equation, forward Euler
                self.u[i,j] = self.u_mdt[i,j] + self.nu*self.dt*(
                                (self.u_mdt[i+1,j]-2.*self.u_mdt[i,j]+self.u_mdt[i-1,j])/self.dx**2
                               +(self.u_mdt[i,j+1]-2.*self.u_mdt[i,j]+self.u_mdt[i,j-1])/self.dy**2)

    #Plot graph
    def plot(self, ax=None):
        if ax == None:
            ax=self.ax

        sm = ax.pcolormesh(self.X,self.Y,self.u.T,vmin=0.,vmax=1.,alpha=0.5)
        draw_grid(ax, nx=self.x.size-1,ny=self.y.size-1, nz=1,
                      xmin=self.xmin,xmax=self.xmax,
                      ymin=self.ymin,ymax=self.ymax)

        #Plot velocity profile offset to the left
        axisloc = self.xmax
        ax.arrow(axisloc,0., 0., 1.,  width=0.0015, color="k", 
                 clip_on=False, head_width=0.012, head_length=0.012)
        ax.arrow(axisloc,0., 1., 0.,  width=0.0015, color="k", 
                 clip_on=False, head_width=0.012, head_length=0.012)
        ax.plot(np.mean(self.u,0)+axisloc,self.y,'g-')
        #ax.set_xlim((0.,2.))

        if self.first_time:
            plt.colorbar(sm)
            self.first_time=False

        plt.pause(0.01)
        ax.cla()

if __name__ == "__main__":

    t0 = 0; tf = 30.; Nsteps = 10000
    time = np.linspace(t0,tf,Nsteps)
    dt = np.mean(np.diff(time))
    uwall = 1.
    
    cfd = CFD(nu=0.575, dt=dt, fig=fig,
              xsize = ncx, ysize = ncy+2,
              xmin = 0., xmax = xl_cfd,
              ymin = 0., ymax = yl_cfd)

    for n,t in enumerate(time):

        print("CFD time = ", n,t)

        #===============================================
        # Call to CPL-LIBRARY goes here to
        # recieve u_MD to set bottom boundary
        #===============================================
        #umd = cpl.recv(u)
        #bottomwall = np.mean(umd)

        #Update CFD
        cfd.set_bc(topwall=uwall, bottomwall=0.)
        cfd.update_time()
        cfd.plot()

        #===============================================
        # Call to CPL-LIBRARY goes here to
        # send u_CFD in constraint region
        #===============================================
        #ucnst = cfd.u[:,7]


