from cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits
from mpi4py import MPI
import numpy as np

# Inheritance for _init() instead of call-backing allows the user
# to define extra useful attributes. Pack and Unpack actions seems more
# naturally to be just callbacks since they keep no state, just build
# the recv_ and sen_ buffer.
class CPLSocket(object) :
    def __init__(self, context=None):
        self.magnitudes = {}
        self._initialised = False
        self._lib = CPL()
        self._realm_comm = None
        self.context = context
        self._realm = CPL.NULL_REALM
        self.nsteps = 0
        # Send and receive internal buffers
        self.recv_buff = []
        self.send_buff = []
        # Flag to check if buffer is ready to be sent
        self.packed = False
        self.packers = {}
        self.unpackers = {}
        self.activePacker = ""
        # Init domain and procs attributes
        self._initialise(0.0, [0,0,0], [0, 0, 0], None)
        self._createComms()

    # _initialise take care of the initialization of the attributes
    # defined in the current class needed to call CPL initialization
    # functions.
    def _initialise(self, dt, npxyz, xyzL, topology): 
        # Domain dimensions
        self.lx = xyzL[0]
        self.ly = xyzL[1]
        self.lz = xyzL[2]
        self.dt = dt 
        self.topology = topology


    def _createComms(self):
        if (not self._initialised) and (self._realm_comm is None):
            self._realm_comm = self._lib.create_comm(self._realm)
        else:
            print "Already created or initialised"
    
    def registerPacker(self, pack_name, pack_callb):
        self.packers[pack_name] = pack_callb


    def registerUnpacker(self, unpack_name, unpack_callb):
        self.unpackers[unpack_name] = unpack_callb


    
    # Init methods
    def init(self, **args):
        if len(args) > 0:
            self._initialise(**args)
        else:
            if self.context is not None:
                self._init()
            else:
                print "Exception there is no context"
                exit()

    # Send methods    

    def send(self, buff=None, pack_name=None):
        if self.context is not None: 
            self._pack(pack_name)
            self.active_packer = pack_name
        else: 
            print "Exception"
            exit()
        self._send(self.send_buff)
        
    
    # Receive methods    
    def receive(self, data):
        return self._receive()

    def receiveDomainData(self, mag):
        if self.context is not None:
            self._unpack(mag)
        else:
            print "Exception"
            exit()
        self._receive()

    # Routines that calls the appropiate pack/unpack callback
    # Their pourpose is to prepare send_buff/recv_buff to be
    # sent/received. 
    def _pack(self, pack_name):
        self.packers[pack_name](self)

    def _unpack(self, unpack_name):
        self.unpackers[unpack_name](self)


    # Abstract methods
    def _init(self):
        pass

    def _send(self, data): 
        pass

    def _receive(self):
        pass
    
 

class ProcTopology(object):
    def __init__(self):
        self._topo_comm = None
        self._topo_coords = None
        self.created = False
        self.nprocs = 0

    def setNumProcs():
        pass

    def createTopo(self, realm_comm):
        pass

    def gridDims():
        pass

    # Check consistency of the number of processes
    def checkProcs(self):
        if (self.topo_comm.Get_size() != self.nprocs):
            print "Non-coherent number of processes"
            exit()

class CartProcTopology(ProcTopology):
    def __init__(self, npx, npy, npz ):
        super(ProcTopology, self).__init__()
        self.npx = npx
        self.npy = npy
        self.npz = npz
        self.setNumProcs()
        self.checkProcs()


    def setNumProcs(self):
        self.nprocs = self.npx * self.npy * self.npz

    def gridDims(self):
        return [self.npx, self.npy, self.npz]

    def createTopo(self, realm_comm, icoords=None):
        if self.topo_comm is None:
            self.topo_comm = realm_comm.Create_cart([self.npx, 
                                                     self.npy, self.npz])
            if icoords is None:
                self._createCartTopo()
            else:
                self._createCartTopoFromCoordsicoords()
            self.created = True

    def _createCartTopo(self):
        for rank in xrange(self.nprocs):
            cart_coords_l = self.topo_comm.Get_coords(rank)
            self.topo_coords[:, rank] = cart_coords_l 
        self.topo_coords += 1
        self.topo_coords = toCPLArray(self.topo_coords, np.int32)
    
    def _createCartTopoFromCoords(self):
        pass


class Packer(object):
    def __init__(self, region, callback):
        self.region = toCPLArray(region, np.int32)
        self.callback = callback

    def __call__(self, socket):
        self.callback(self.region, socket)


class CPLGrid(object):
    def __init_(self, topo_comm, xg=[], yg=[], zg=[]):
        self.created = False
        if (xg and yg and zg):
            self.xg = toCPLArray(xg, np.float64)
            self.yg = toCPLArray(yg, np.float64)
            self.zg = toCPLArray(yg, np.float64)
            self.created = True
        self.topo_comm = topo_comm
        self.proc_bounds = []
        self.glob_bounds = []

    def createGrid():
        pass

    def getCPLProcBounds(rank=None):
        pass

    def computeCPLProcsBounds(self):
        nprocs = self.topo_comm.Get_size()
        pb = self.procs_bounds
        pb = np.zeros((6, nprocs), order='F', dtype=np.int32)
        for rank in xrange(nprocs):
            pb[:, rank] = self.getCPLProcBounds(rank)
            

    def computeCPLGlobBounds(self):
        pass

class CPLCartGrid(CPLGrid):
    def __init_(self, topo_comm, xg=None, yg=None, zg=None, ncxyz=None, lxyz=None):
        super(CPLCartGrid, self).__init__(topo_comm, xg, yg, zg)
        if (ncxyz is not None) and (lxyz is not None):
            self.ncx = ncxyz[0]
            self.ncy = ncxyz[1]
            self.ncz = ncxyz[2]

            if not self.created:
                self.createGrid(ncxyz, lxyz)
        else:
            self.created = False
            

    def createGrid(self, ncxyz, lxyz):
        if ncxyz and lxyz:
            self.ncx = ncxyz[0]
            self.ncy = ncxyz[1]
            self.ncz = ncxyz[2]
            dx = lxyz[0]/self.ncx
            dy = lxyz[1]/self.ncy
            dz = lxyz[2]/self.ncz
            [self.ncxl, self.ncyl, self.nczl] = ncxyz / lxyz
            [self.xg, self.yg, self.zg] = create_CPL_cart_3Dgrid(self.ncx, 
                                                                self.ncy, 
                                                                self.ncz, 
                                                                dx, dy, dz)
            self.computeGlobBounds()
            self.computeCPLProcsBounds()
            self.created = True
        else:
            print "Exception not enough params"

    def getCPLProcBounds(self, rank):
	cart_coords = self.topo_comm.Get_coords(rank)
        [x, y, z] = cart_coords
        iTmin = x*self.ncxl + 1
	iTmax = iTmin + self.ncxl - 1
        jTmin = y*self.ncyl + 1
	jTmax = jTmin + self.ncyl - 1
        kTmin = z*self.nczl + 1
	kTmax = kTmin + self.nczl - 1
        return [iTmin, iTmax, jTmin, jTmax, kTmin, kTmax]

    def computeCPLGlobBounds(self):
        self.glob_bounds = np.array([[1, self.ncx], [1, self.ncy], [1, self.ncz]],
                                    order='F', dtype=np.int32)


            


class CPLSocketMD(CPLSocket):
    def __init__(self, context=None):
        super(CPLSocketMD, self).__init__(context)
        self.initial_step = 0
        self._realm = CPL.MD_REALM

    def _initialise(self, dt, npxyz, xyzL, topology): 
        super(CPLSocketMD, self)._initialise(dt, npxyz, xyzL, topology)


    def _send(self, data):
        # This should not happen since packed() is called in send() but
        # deffensive programming is never a bad thing.
        if self.packed:
            olap_limits = get_olap_limits(self._lib)
#            self._lib.scatter(self.recv_buff, olap_limits, self.send_buff)
            self.packed = False
        else:
            print "Exception: pack() has not been called"


    def _receive(self, data):
        olap_limits = get_olap_limits(self._lib)
#        self._lib.gather(self.send_buff, olap_limits, self.recv_buff)

    def _init(self):
        self._initMD()
        npxyz = self.topology.gridDims()
        xyzL = [self.lx, self.ly, self.lz]
        dummy_density = 1.0
        self.nsteps, self.initialstep = self._lib.md_init(self.dt, 
                                                self.topology.topo_comm, 
                                                self.topology.topo_coords, 
                                                npxyz, xyzL, dummy_density)
    def _initMD(self):
        pass

class CPLSocketCFD(CPLSocket):
    def __init__(self, context=None):
        super(CPLSocketCFD, self).__init__(context)
        self.ncx = 0
        self.ncy = 0
        self.ncz = 0
        self._realm = CPL.CFD_REALM


    def _initialise(self, dt, npxyz, xyzL, topology, grid): 
        super(CPLSocketCFD, self)._initialise(dt, npxyz, xyzL, topology)
        self.grid = grid


    def _send(self, data):
        if self.packed:
            limits = self.packers[self.active_packer].region
#           self._lib.scatter(self.recv_buff, limits, self.send_buff)
            self.packed = False
        else:
            print "Exception: pack() has not been called"

    def _receive(self, data):
        limits = self.unpackers[self.active_packer].region
#        self._lib.gather(self.send_buff, limits, self.recv_buff)

    def _init(self):
        self._initCFD()
        npxyz = [self.topology.npx, self.topology.npy, self.topology.npz]
        xyzL = [self.lx, self.ly, self.lz]
        ncxyz = [self.grid.ncx, self.grid.ncy, self.grid.ncz]
        dummy_density = 1.0
        pb = self.proc_bounds
        [iTmin, iTmax] = pb[0:2, :]
        [jTmin, jTmax] = pb[2:4, :]
        [kTmin, kTmax] = pb[4:6, :]
        ijkcmin = self.glob_bounds[:, 0]
        ijkcmax = self.glob_bounds[:, 1]
        self._lib.cfd_init(self.dt, self.topology.topo_comm, 
                            self.topology.topo_coords, npxyz, xyzL, ncxyz,
                            dummy_density, ijkcmax, ijkcmin, iTmin, iTmax, 
                            jTmin, jTmax, kTmin,kTmax, self.grid.xg, 
                            self.grid.yg, self.grid.zg)



def toCPLArray(arr, arr_type):
    if type(arr_type) != np.ndarray:
        return np.array(arr, order='F', dtype=arr_type)



