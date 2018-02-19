#import numpy.ctypeslib as ctl
#import ctypes
from overlap import *
import matplotlib.pyplot as plt
import numpy as np

#Load c++ interp library
#libname = 'interplib.so'
#libdir = './'
#lib=ctl.load_library(libname, libdir)
#sphere_cube_overlap = lib.sphere_cube_overlap
#sphere_cube_overlap.argtypes = [ctypes.c_double]*10
#sphere_cube_overlap.restype = ctypes.c_double

#Define some values
xi = 0.0; yi=0.5; zi=0.0; rp = 1.;
xb = -1.0; yb=-1.0; zb=-1.0;
xt =  1.0; yt= 1.0; zt= 1.0;
r = [xi, yi, zi]
box = [xb, yb, zb, xt, yt, zt]

#Control volume function
xmin = min(xi-rp, xb)
xmax = max(xi+rp, xt)
ymin = min(yi-rp, yb)
ymax = max(yi+rp, yt)
zmin = min(zi-rp, zb)
zmax = max(zi+rp, zt)
xyzminmax = [xmin, ymin, zmin, xmax, ymax, zmax]

Nx = 200
Ny = 200
Nz = 200
x, y, z = np.meshgrid(np.linspace(xmin, xmax, Nx), 
                      np.linspace(ymin, ymax, Ny),
                      np.linspace(zmin, zmax, Nz), indexing='ij')

CVsphr = CV(x, y, z, xb, xt, yb, yt, zb, zt)*sphere(xi, yi, zi, rp, x, y, z)

#Get sum of area
dx =(xmax-xmin)/Nx 
dy =(ymax-ymin)/Ny 
dz =(zmax-zmin)/Nz 
dV = dx*dy*dz
xloc = int((xi-xmin)/dx)
yloc = int((yi-ymin)/dy)
zloc = int((zi-zmin)/dz)

print(xyzminmax, box, r, rp, xloc, yloc, zloc)
sco = get_overlap(r, rp, box)
print("Sum of area ", (4./3.)*np.pi*rp**3, np.sum(CVsphr)*dV, sco,  sco/((4./3.)*np.pi*rp**3))
fig, ax = plt.subplots(3,1)

ax[0].pcolormesh(x[:,:,zloc], y[:,:,zloc], CVsphr[:,:,zloc])
drawbox(ax[0], xb, xt, yb, yt)
ax[1].pcolormesh(y[xloc,:,:], z[xloc,:,:], CVsphr[xloc,:,:])
drawbox(ax[1], yb, yt, zb, zt)
ax[2].pcolormesh(x[:,yloc,:], z[:,yloc,:], CVsphr[:,yloc,:])
drawbox(ax[2], xb, xt, zb, zt)

ax[0].plot(xi, yi, 'gx')
ax[1].plot(yi, zi, 'gx')
ax[2].plot(xi, zi, 'gx')

plt.show()



