import numpy as np
import matplotlib.pyplot as plt

def draw_lines(ax, x, y, lc, lw=1., XKCD_plots=False):

    if XKCD_plots:
        from XKCD_plot_generator import xkcd_line

    for n in range(x.size):
        if XKCD_plots:
            xp, yp = xkcd_line([x[n], x[n]],  [y[0], y[-1]])
        else:
            xp, yp = [x[n], x[n]],  [y[0], y[-1]]      
        ax.plot(xp, yp, '-'+lc,lw=lw,alpha=0.5)
    for n in range(y.size):
        if XKCD_plots:
            xp, yp = xkcd_line([x[0], x[-1]], [y[n], y[n]])
        else:
            xp, yp = [x[0], x[-1]], [y[n], y[n]]
        ax.plot(xp, yp, '-'+lc,lw=lw,alpha=0.5)


def draw_grid(ax, nx, ny, nz,
              px=1, py=1, pz=1,
              xmin=0.,ymin=0.,zmin=0.,
              xmax=1.,ymax=1.,zmax=1.,
              ind=0., fc='r', lc='k',
              offset=0., label="",
              nodes=False,
              XKCD_plots=False, **kwargs):

    ymin += offset
    ymax += offset

    x = np.linspace(xmin, xmax, nx+1)
    y = np.linspace(ymin, ymax, ny+1)
    z = np.linspace(zmin, zmax, nz+1)
    if nodes:
        X,Y,Z = np.meshgrid(x,y,z)
        self.ax.plot(X[:,:,ind],Y[:,:,ind],'s'+fc,alpha=0.5)

    draw_lines(ax, x, y, lc=lc,XKCD_plots=XKCD_plots)

    xp = np.linspace(xmin, xmax, px+1)
    yp = np.linspace(ymin, ymax, py+1)
    zp = np.linspace(zmin, zmax, pz+1)
    if nodes:
        X,Y,Z = np.meshgrid(xp,yp,zp)
        ax.plot(X[:,:,ind],Y[:,:,ind],'s'+fc,alpha=0.5)

    draw_lines(ax, xp, yp, lc='k',lw=2,XKCD_plots=XKCD_plots)

    #Annotate the processors
    dxp = np.mean(np.diff(xp))
    dyp = np.mean(np.diff(yp))
    #Make some effort to ensure the text fits on a normal screen
    if px > 8 or py > 8:
        fontsize = 8
    else:
        fontsize = 12
    for i,px in enumerate(xp[1:]):
        for j,py in enumerate(yp[1:]):
            ax.text(px-.5*dxp,py-.5*dyp,
                    label+"$ Proc [" + str(i) + "," + str(j) + "]$",
                    horizontalalignment="center",color=lc,
                    fontsize=fontsize)

