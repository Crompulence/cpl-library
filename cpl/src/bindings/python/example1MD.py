from cplsocket import *

# Custom socket to give an example of how to handle self.context
class MyCartSocketMD(CartSocketMD):
    def __init__(self, context):
        CartSocketMD.__init__(self, context)
        # Do stuff with myvar1 and myvar2, define new
        # attributes from the self.context, etc. 

    def _initMD(self):
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
            "npx" : 3, "npy" : 3, "npz" : 3,
            "dt" : 0.1}


md_socket = MyCartSocketMD(MDcontext)
md_socket.init()




