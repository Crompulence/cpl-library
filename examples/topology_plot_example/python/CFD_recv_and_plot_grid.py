import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from cplpy import CPL
from draw_grid import draw_grid

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
nprocs_realm = CFD_COMM.Get_size()

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
NProcs = np.product(npxyz)
xyzL = np.array([10.0, 10.0, 10.0], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([16, 6, 16], order='F', dtype=np.int32)

if (nprocs_realm != NProcs):
    print(("Non-coherent number of processes in CFD ", nprocs_realm,
            " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2]))
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

# recv data to plot
olap_limits = CPL.get_olap_limits()
portion = CPL.my_proc_portion(olap_limits)
[ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
recv_array = np.zeros((1, ncxl, ncyl, nczl), order='F', dtype=np.float64)
recv_array, ierr = CPL.recv(recv_array, olap_limits)

#Plot output
fig, ax = plt.subplots(1,1)

# === Plot both grids ===
dx = CPL.get("xl_cfd")/float(CPL.get("ncx"))
dy = CPL.get("yl_cfd")/float(CPL.get("ncy"))
dz = CPL.get("zl_cfd")/float(CPL.get("ncz"))
ioverlap = (CPL.get("icmax_olap")-CPL.get("icmin_olap")+1)
joverlap = (CPL.get("jcmax_olap")-CPL.get("jcmin_olap")+1)
koverlap = (CPL.get("kcmax_olap")-CPL.get("kcmin_olap")+1)
xoverlap = ioverlap*dx
yoverlap = joverlap*dy
zoverlap = koverlap*dz

print((CPL.get("xl_cfd")))

#Plot CFD and coupler Grid
draw_grid(ax, 
          nx=CPL.get("ncx"),
          ny=CPL.get("ncy"),
          nz=CPL.get("ncz"),
          px=CPL.get("npx_cfd"),
          py=CPL.get("npy_cfd"),
          pz=CPL.get("npz_cfd"),
          xmin=CPL.get("x_orig_cfd"),
          ymin=CPL.get("y_orig_cfd"),
          zmin=CPL.get("z_orig_cfd"),
          xmax=(CPL.get("icmax_olap")+1)*dx,
          ymax=CPL.get("yl_cfd"),
          zmax=(CPL.get("kcmax_olap")+1)*dz,
          lc = 'r',
          label='CFD')

#Plot MD domain
draw_grid(ax, nx=1, ny=1, nz=1,
          px=CPL.get("npx_md"),
          py=CPL.get("npy_md"),
          pz=CPL.get("npz_md"),
          xmin=CPL.get("x_orig_md"),
          ymin=-CPL.get("yl_md")+yoverlap,
          zmin=CPL.get("z_orig_md"),
          xmax=(CPL.get("icmax_olap")+1)*dx,
          ymax=yoverlap,
          zmax=(CPL.get("kcmax_olap")+1)*dz,
          label='MD')

#Plot some random molecules
#ax.plot(np.random.random(100)*(CPL.get("xl_md")),
#           np.random.random(100)*(CPL.get("yl_md"))-CPL.get("yl_md")+yoverlap,
#           'ob',alpha=0.5)

#Plot x component on grid
x = np.linspace(CPL.get("x_orig_cfd")+.5*dx,xoverlap-.5*dx,ioverlap)
z = np.linspace(CPL.get("z_orig_cfd")+.5*dz,zoverlap-.5*dz,koverlap)
for j in range(joverlap):
    ax.plot(x, 0.5*dy*(recv_array[0,:,j,0]+1.+2*j), 's-')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.show()

CPL.finalize()
MPI.Finalize()




