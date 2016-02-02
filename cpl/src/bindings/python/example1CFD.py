from cplsocket import *

# Custom socket to give an example of how to handle self.context
class MyCartSocketCFD(CartSocketCFD):
    def __init__(self, context, myvar1, myvar2):
        CartSocketCFD.__init__(self, context)
                # Do stuff with myvar1 and myvar2, define new
                # attributes from the self.context, etc. 

    def _initCFD(self):
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
            "ncx" : 18, "ncy" : 18, "ncz" : 18,
            "npx" : 3, "npy" : 3, "npz" : 3,
            "dt" : 0.1}

var1 = "var1"
var2 = "var2"
cfd_socket = MyCartSocketCFD(CFDcontext, var1, var2)
cfd_socket.init()

# Packing velocities
def pack_vel_cb(:




