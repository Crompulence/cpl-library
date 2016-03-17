from cplpy.cplsocket import CartSocketMD, CartProcTopology
from cplpy.cplpacker import Packer, Unpacker
from cplpy.cpl import CPL
import numpy as np
import matplotlib.pyplot as plt

from draw_grid import draw_grid

# Custom socket to give an example of how to handle self.context
class MyCartSocketMD(CartSocketMD):
    def __init__(self, context):
        CartSocketMD.__init__(self, context)

    def _initMD(self):
        self.nsteps = self.context["nsteps"]
        self.lx = self.context["lx"]	
        self.ly = self.context["ly"]	
        self.lz = self.context["lz"]	
        self.dt = self.context["dt"]
        npx = self.context["npx"]
        npy = self.context["npy"]
        npz = self.context["npz"]
        self.topology = CartProcTopology(npx, npy, npz)
        self.topology.createTopo(self._realm_comm)


# Implement and register packing and unpacking functions
MDcontext = {"lx" : 10, "ly" : 5, "lz" : 10,
            "npx" : 4, "npy" : 2, "npz" : 1,
            "dt" : 0.1, "nsteps":100}

root = 0
md_socket = MyCartSocketMD(MDcontext)
md_socket.init()
lib = md_socket._lib

#Get all info from cpl-library
inoverlap = lib.get('overlap')
ncx = lib.get("ncx")
ncy = lib.get("ncy")
ncz = lib.get("ncz")
xl_md = lib.get("xl_md")
yl_md = lib.get("yl_md")
zl_md = lib.get("zl_md")
xl_cfd = lib.get("xl_cfd")
yl_cfd = lib.get("yl_cfd")
zl_cfd = lib.get("zl_cfd")
npx_md = lib.get("npx_md")
npy_md = lib.get("npy_md")
npz_md = lib.get("npz_md")
npx_cfd = lib.get("npx_cfd")
npy_cfd = lib.get("npy_cfd")
npz_cfd = lib.get("npz_cfd")
xmin_olap = lib.get("icmin_olap")
xmax_olap = lib.get("icmax_olap")
ymin_olap = lib.get("jcmin_olap")
ymax_olap = lib.get("jcmax_olap")
zmin_olap = lib.get("kcmin_olap")
zmax_olap = lib.get("kcmax_olap")
dx = xl_cfd/float(ncx)
dy = yl_cfd/float(ncy)
dz = zl_cfd/float(ncz)
xmin_olap_md = (xmin_olap-1)*dx
ymin_olap_md = (ymin_olap-1)*dy
zmin_olap_md = (zmin_olap-1)*dz
xmax_olap_md = xmax_olap*dx
ymax_olap_md = ymax_olap*dy
zmax_olap_md = zmax_olap*dz
xoverlap = xmax_olap_md-xmin_olap_md #(xmax_olap-xmin_olap+1)*dx
yoverlap = ymax_olap_md-ymin_olap_md #(ymax_olap-ymin_olap+1)*dy
zoverlap = zmax_olap_md-zmin_olap_md #(zmax_olap-zmin_olap+1)*dz

# === Plot both grids ===
#Get overlap region
comm = md_socket.topology._topo_comm
rank = comm.Get_rank()
COMM_size = comm.Get_size()
pcoords=comm.Get_coords(rank)
if rank == root:

    #Plot
    fig, ax = plt.subplots(1,1)
    draw_grid(ax, nx=ncx,
                  ny=ncy,
                  nz=ncz,
                  px=npx_cfd,
                  py=npy_cfd,
                  pz=npz_cfd,
                  xmin=0.,ymin=0.,zmin=0.,
                  xmax=xl_cfd,
                  ymax=yl_cfd,
                  zmax=zl_cfd,
                  lc = 'r',
                  label='CFD')

    draw_grid(ax, nx=1, ny=1, nz=1,
                  px=npx_md,
                  py=npy_md,
                  pz=npz_md,
                  xmin=-xl_md+xoverlap,
                  ymin=-yl_md+yoverlap,
                  zmin=-zl_md+zoverlap,
                  xmax=xoverlap,
                  ymax=yoverlap,
                  zmax=zoverlap,
                  label='MD')

    #Plot some random molecules
    ax.plot(np.random.random(100)*(xl_md),
            np.random.random(100)*(yl_md)-yl_md+yoverlap,
            'ob',alpha=0.5)

    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')



#Recv data
icmin = 1; icmax = ncx
jcmin = 1; jcmax = 4
kcmin = 1; kcmax = ncz
limits = np.array([icmin, icmax, jcmin, jcmax, kcmin, kcmax],dtype=np.int32)

A = np.zeros([4,icmax-icmin+1,jcmax-jcmin+1,kcmax-kcmin+1])
A, ierr = lib.recv(A, *limits)

#Plot recieved data
x = np.linspace(.5*dx,xl_cfd-.5*dx,ncx)
z = np.linspace(.5*dz,zl_cfd-.5*dz,ncz)
Agather = comm.gather(A, root=root)
inoverlapgather = comm.gather(inoverlap, root=root)
colors = [["g","c"],["c","g"]]
if rank == root:
    nrank = -1
    for ib in range(npx_md):
        for jb in range(npy_md):
            portion = lib.proc_portion(np.array([ib+1,jb+1,1],dtype=np.int32), lib.MD_REALM, limits)
            #extents = lib.proc_extents(np.array([ib+1,jb+1,1],dtype=np.int32), lib.MD_REALM)
            nxp = portion[1]-portion[0]+1
            nyp = portion[3]-portion[2]+1
            nrank += 1
            print("MD shape", nrank, inoverlapgather[nrank], portion, Agather[nrank].shape, np.sum(Agather[nrank]))
            if inoverlapgather[nrank]:
                for j in range(portion[2]-1,portion[3]):
                    print(nrank,j,j-portion[2])
                    ax.plot(x[portion[0]-1:portion[1]], j*dy+0.5*dy*(Agather[nrank][0,0:nxp,j-portion[2]+1,0]+1.),'-s',alpha=0.5, color=colors[ib%2][jb%2])

plt.show()
