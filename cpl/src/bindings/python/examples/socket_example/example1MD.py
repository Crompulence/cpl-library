from cplpy.cplsocket import CartSocketMD, CartProcTopology
from cplpy.cplpacker import Packer, Unpacker
from cplpy.cpl import CPL
import numpy as np

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

MDcontext = {"lx" : 10, "ly" : 10, "lz" : 10,
            "npx" : 3, "npy" : 3, "npz" : 1,
            "dt" : 0.1, "nsteps":100}


md_socket = MyCartSocketMD(MDcontext)
md_socket.init()
#md_socket.comm_mode = CPL.GATHER_SCATTER
md_socket.comm_mode = CPL.SEND_RECEIVE

# Packing velocities
def pack_vel_cb(socket):
    socket.send_buff *= -socket.topology._topo_comm.Get_rank()

vel_reg = md_socket.olap_region
pack_vel = Packer("velocity", vel_reg, pack_vel_cb, 3)
md_socket.registerPacker(pack_vel)

# Packing stress
def unpack_stress_cb(socket):
    if np.any(socket.recv_buff > -1):
        print socket.recv_buff

stress_reg = md_socket.olap_region
unpack_stress = Unpacker("stress", stress_reg, unpack_stress_cb, 9)
md_socket.registerUnpacker(unpack_stress)

# Main Loop
md_socket.receive("stress")
md_socket.send("velocity")
