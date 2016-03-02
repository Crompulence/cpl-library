from cplpy.cplsocket import CartSocketCFD, CartProcTopology, CartGrid
from cplpy.cplpacker import Packer, Unpacker
from cplpy.cpl import CPL
import numpy as np
import matplotlib.pyplot as plt

from draw_grid import draw_grid




# Custom socket to give an example of how to handle self.context
class MyCartSocketCFD(CartSocketCFD):
    def __init__(self, context):
        CartSocketCFD.__init__(self, context)

    def _initCFD(self):
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
        ncx = self.context["ncx"]
        ncy = self.context["ncy"]
        ncz = self.context["ncz"]
        lx = self.context["lx"]
        ly = self.context["ly"]
        lz = self.context["lz"]
        self.grid = CartGrid([ncx, ncy, ncz], [lx, ly, lz])


lib = CPL()

CFDcontext = {"lx" : 10, "ly" : 10, "lz" : 10,
              "ncx" : 32, "ncy" : 32, "ncz" : 1,
              "npx" : 1, "npy" : 1, "npz" : 1,
              "dt" : 0.2, "nsteps":100}

cfd_socket = MyCartSocketCFD(CFDcontext)
cfd_socket.init()
#cfd_socket.comm_mode = CPL.GATHER_SCATTER
cfd_socket.comm_mode = CPL.SEND_RECEIVE

# Unpacking velocities
def unpack_vel_cb(socket):
    if np.any(socket.recv_buff > -1):
        print socket.recv_buff

vel_reg = cfd_socket.olap_region
unpack_vel = Unpacker("velocity", vel_reg, unpack_vel_cb, 3)
cfd_socket.registerUnpacker(unpack_vel)

# Packing stress
def pack_stress_cb(socket):
    socket.send_buff *= -socket.topology._topo_comm.Get_rank()

stress_reg = cfd_socket.olap_region
pack_stress = Packer("stress", stress_reg, pack_stress_cb, 9)
cfd_socket.registerPacker(pack_stress)


# === Plot both grids ===

#Get overlap region
if cfd_socket.topology._topo_comm.Get_rank() == 0:

    #Get all info from cpl-library
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
    dx = xl_cfd/float(CFDcontext['ncx'])
    dy = yl_cfd/float(CFDcontext['ncy'])
    dz = zl_cfd/float(CFDcontext['ncz'])
    xmin_olap_md = (xmin_olap-1)*dx
    ymin_olap_md = (ymin_olap-1)*dy
    zmin_olap_md = (zmin_olap-1)*dz
    xmax_olap_md = xmax_olap*dx
    ymax_olap_md = ymax_olap*dy
    zmax_olap_md = zmax_olap*dz
    xoverlap = xmax_olap_md-xmin_olap_md#(xmax_olap-xmin_olap+1)*dx
    yoverlap = ymax_olap_md-ymin_olap_md#(ymax_olap-ymin_olap+1)*dy
    zoverlap = zmax_olap_md-zmin_olap_md#(zmax_olap-zmin_olap+1)*dz

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

                  

# Main Loop
cfd_socket.send("stress")
cfd_socket.receive("velocity")

print(cfd_socket.recv_buff)
#ax.plot(cfd_socket.recv_buff)
plt.show()
