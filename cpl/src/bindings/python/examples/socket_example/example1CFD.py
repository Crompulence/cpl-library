from cplpy.cplsocket import CartSocketCFD, CartProcTopology, CartGrid
from cplpy.cplpacker import Packer, Unpacker
from cplpy.cpl import CPL
import numpy as np

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



CFDcontext = {"lx" : 10, "ly" : 10, "lz" : 10,
            "ncx" : 6, "ncy" : 6, "ncz" : 1,
            "npx" : 3, "npy" : 3, "npz" : 1,
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

# Main Loop
cfd_socket.send("stress")
cfd_socket.receive("velocity")
