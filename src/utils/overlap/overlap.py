import numpy.ctypeslib as ctl
import ctypes
import os 
import numpy as np
import matplotlib.pyplot as plt

def H(x):
    return 1 * (x > 0)

def sphere(xi, yi, zi, rp, x, y, z):
    return H(rp - np.sqrt(((x-xi)**2 + (y-yi)**2+ (z-zi)**2)) )

def CV(x, y, z, xb, xt, yb, yt, zb, zt):
    return ( (H(xt - x) - H(xb - x))
            *(H(yt - y) - H(yb - y))
            *(H(zt - z) - H(zb - z)))


def drawbox(ax, xp=0.0, xm=-1.0, yp=0.9, ym=-0.1):

    ax.plot((xm,xp), (ym, ym), 'k-', lw=3.)
    ax.plot((xp,xp), (ym, yp), 'k-', lw=3.)
    ax.plot((xm,xp), (yp, yp), 'k-', lw=3.)
    ax.plot((xm,xm), (ym, yp), 'k-', lw=3.)


def plot_porosity(ax, box, r, rp, dxyz, Nx = 100, Ny = 100, Nz = 100):

    xb = box[0]; yb = box[1]; zb = box[2]
    xt = box[3]; yt = box[4]; zt = box[5]
    drawbox(ax[0], xb, xt, yb, yt)
    drawbox(ax[1], yb, yt, zb, zt)
    drawbox(ax[2], xb, xt, zb, zt)

    xi, yi, zi = r
    dx, dy, dz = dxyz
    dV = dx*dy*dz

    xmin = min(xi-rp, xb)
    xmax = max(xi+rp, xt)
    ymin = min(yi-rp, yb)
    ymax = max(yi+rp, yt)
    zmin = min(zi-rp, zb)
    zmax = max(zi+rp, zt)

    xloc = int((xi-xmin)/dx)
    yloc = int((yi-ymin)/dy)
    zloc = int((zi-zmin)/dz)

    dxm =(xt-xb)/Nx 
    dym =(yt-yb)/Ny 
    dzm =(zt-zb)/Nz 
    dVm = dxm*dym*dzm

    x, y, z = np.meshgrid(np.linspace(xb, xt, Nx), 
                          np.linspace(yb, yt, Ny),
                          np.linspace(zb, zt, Nz), indexing='ij')

    #Control volume function
    CVsphr = CV(x, y, z, xb, xt, yb, yt, zb, zt)*sphere(xi, yi, zi, rp, x, y, z)

    ax[0].pcolormesh(x[:,:,zloc], y[:,:,zloc], CVsphr[:,:,zloc], alpha=0.2)
    ax[1].pcolormesh(y[xloc,:,:], z[xloc,:,:], CVsphr[xloc,:,:], alpha=0.2)
    ax[2].pcolormesh(x[:,yloc,:], z[:,yloc,:], CVsphr[:,yloc,:], alpha=0.2)

    #print("Cell", i,j,k,np.sum(CVsphr)*dVm,cell[1+i,1+j,1+k], np.abs((np.sum(CVsphr)*dVm-cell[1+i,1+j,1+k])/cell[1+i,1+j,1+k]))
    #assert np.abs((np.sum(CVsphr)*dVm-cell[1+i,1+j,1+k])/cell[1+i,1+j,1+k]) < 0.2

#Load c++ interp library
def get_overlap(r, rp, box):
    """
        Get fraction of overlap in cell by calling c++ code
        r -- Position of particle
        rp -- radius of particle
        box -- xbottom, ybottom, zbottom, xtop, ytop, ztop
    """

    libname = 'overlaplib.so'
    libdir = os.path.dirname(os.path.realpath(__file__)) #'./'
    lib=ctl.load_library(libname, libdir)
    sphere_cube_overlap = lib.sphere_cube_overlap
    Nargs = len(r)+1+len(box)
    sphere_cube_overlap.argtypes = [ctypes.c_double]*Nargs
    sphere_cube_overlap.restype = ctypes.c_double

    return sphere_cube_overlap(r[0], r[1], r[2], rp, 
                               box[0], box[1], 
                               box[2], box[3], 
                               box[4], box[5])


def get_neighbour_cells(cellcentre, dxyz, r, rp, plot=False):

    """
        Assumes particle in central cell and we take the
        27 cells around current cell.
    """

    maxside = 1
    for side in dxyz:
        maxside = int(max(np.ceil(rp/side), maxside))
        print(side, np.ceil(rp/side), rp, maxside)

    cell = np.zeros([1+2*maxside,1+2*maxside,1+2*maxside])

    #Map particle position so cell centre is zero
    rm = r - cellcentre

    if plot:
        fig, ax = plt.subplots(3,1)

    for i in range(-maxside,maxside+1):
        for j in range(-maxside,maxside+1):
            for k in range(-maxside,maxside+1):
                box[0] = i*dxyz[0]-0.5*dxyz[0]
                box[1] = j*dxyz[1]-0.5*dxyz[1]
                box[2] = k*dxyz[2]-0.5*dxyz[2]
                box[3] = i*dxyz[0]+0.5*dxyz[0]
                box[4] = j*dxyz[1]+0.5*dxyz[1]
                box[5] = k*dxyz[2]+0.5*dxyz[2]

                cell[maxside+i,maxside+j,maxside+k] = get_overlap(rm, rp, box)
                if cell[maxside+i,maxside+j,maxside+k] > 1e-5:
                    if plot:
                        plot_porosity(ax, box, r, rp, dxyz)

    if plot:
        plt.show()

    return cell


def get_27_neighbour_cells(cellcentre, dxyz, r, rp, plot=False):

    """
        Assumes particle in central cell and we take the
        27 cells around current cell.
    """

    cell = np.zeros([3,3,3])

    #Map particle position so cell centre is zero
    rm = r - cellcentre

    if plot:
        fig, ax = plt.subplots(3,1)

    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                box[0] = i*dxyz[0]-0.5*dxyz[0]
                box[1] = j*dxyz[1]-0.5*dxyz[1]
                box[2] = k*dxyz[2]-0.5*dxyz[2]
                box[3] = i*dxyz[0]+0.5*dxyz[0]
                box[4] = j*dxyz[1]+0.5*dxyz[1]
                box[5] = k*dxyz[2]+0.5*dxyz[2]

                cell[1+i,1+j,1+k] = get_overlap(rm, rp, box)

                if cell[1+i,1+j,1+k] > 1e-5:
                    if plot:
                        plot_porosity(ax, box, r, rp, dxyz)

    if plot:
        plt.show()

    return cell

    #Get box coordinates
#    for ixyz in range(3):
#        box[ixyz]   = - 0.5*dxyz[ixyz]
#        box[ixyz+3] = + 0.5*dxyz[ixyz]

    #Add overlap to cell
    #cell[1,1,1] = get_overlap(rm, rp, box)

#Define some values
if __name__ == "__main__":
    xi = 0.5; yi=0.5; zi=0.5; rp = 1.0;
    xb = -1.0; yb=-1.0; zb=-1.0;
    xt =  1.0; yt= 1.0; zt= 1.0;

    libname = 'interplib.so'
    libdir = './'
    lib=ctl.load_library(libname, libdir)
    sphere_cube_overlap = lib.sphere_cube_overlap
    sphere_cube_overlap.argtypes = [ctypes.c_double]*10
    sphere_cube_overlap.restype = ctypes.c_double

    #Get fraction
    result = sphere_cube_overlap(xi, yi, zi, rp, xb, yb, zb, xt, yt, zt)
    r = [xi, yi, zi]
    box = [xb, yb, zb, xt, yt, zt]
    result2 = get_overlap(r, rp, box)
    print(result, result2)  

    #Try random boxes and particles
#    for n in range(10000):
#        xi, yi, zi = np.random.randn(3)
#        xb, yb, zb, dx, dy, dz = np.random.randn(6)
#        xt=xb+np.abs(dx)
#        yt=yb+np.abs(dy)
#        zt=zb+np.abs(dz)
#        result = sphere_cube_overlap(xi, yi, zi, rp, xb, yb, zb, xt, yt, zt)
#        r = [xi, yi, zi]
#        box = [xb, yb, zb, xt, yt, zt]
#        result2 = get_overlap(r, rp, box)
#        print(result, result2)  
#        assert np.abs(result - result2) < 1e-3

    #Test random particles with fixed block of 27 neighbour cells
    xi = 0.0; yi=0.; zi=0.; rp = 20.;
    xb = -1.0; yb=-1.0; zb=-1.0;
    xt =  1.0; yt= 1.0; zt= 1.0;

    r = np.array([xi, yi, zi])
    V = (4./3.)*np.pi*rp**3
    dxyz = np.array([xt-xb, yt-yb, zt-yb])
    Vcell = dxyz[0]*dxyz[1]*dxyz[2]
    cellcentre = np.array([0.5*(xt+xb), 0.5*(yt+yb), 0.5*(zt+yb)])
    eps = get_neighbour_cells(cellcentre=cellcentre, dxyz=dxyz, r=r, rp=rp, plot=False)/V
    print(eps, np.sum(eps))



    from matplotlib.widgets import Slider

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    axzloc = plt.axes([0.2, 0.1, 0.65, 0.03])
    zloc = int(eps.shape[2]/2.)
    szloc = Slider(axzloc, 'zloc', 0, eps.shape[2], valinit=zloc, valfmt='%0.0f')
    cm = ax.pcolormesh(eps[:,:,zloc])

    def update(val):
        zloc = int(szloc.val)
        cm.set_array(eps[:,:,zloc-1].ravel())
        fig.canvas.draw_idle()

    szloc.on_changed(update)
    plt.show()
        

#    for n in range(10000):
#        xi, yi, zi = np.random.randn(3)
#        r = np.array([xi, yi, zi])
#        V = (4./3.)*np.pi*rp**3
#        dxyz = np.array([xt-xb, yt-yb, zt-yb])
#        Vcell = dxyz[0]*dxyz[1]*dxyz[2]
#        cellcentre = np.array([0.5*(xt+xb), 0.5*(yt+yb), 0.5*(zt+yb)])
#        eps = get_neighbour_cells(cellcentre=cellcentre, dxyz=dxyz, r=r, rp=rp, plot=True)/V

#        print(eps, np.sum(eps))


