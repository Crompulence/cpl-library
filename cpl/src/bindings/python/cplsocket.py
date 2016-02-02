from __future__ import division
from cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits, toCPLArray
import numpy as np

# Inheritance for _init() instead of call-backing allows the user
# to define extra useful attributes. Pack and Unpack actions seems more
# naturally to be just callbacks since they keep no state, just build
# the recv_ and sen_ buffer.
class CPLSocket(object) :
    _realm = CPL.NULL_REALM
    def __init__(self, context=None):
        self._initialised = False
        self._lib = CPL()
        self._realm_comm = None
        self.context = context
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
        self.lx = 0
        self.ly = 0
        self.lz = 0
        self.dt = 0.0
        self.topology = None
        self.proc_region = None
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
    
    def registerPacker(self, packer):
        self.packers[packer.name] = packer.callback


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

    def setNumProcs(self):
        pass

    def createTopo(self, realm_comm):
        pass

    def gridDims(self):
        pass

    # Check consistency of the number of processes
    def checkProcs(self):
        if (self._topo_comm.Get_size() != self.nprocs):
            print "Non-coherent number of processes"
            exit()
# Add creation to constructor
class CartProcTopology(ProcTopology):
    def __init__(self, npx, npy, npz):
        ProcTopology.__init__(self)
        self.npx = npx
        self.npy = npy
        self.npz = npz
        self.setNumProcs()
        self._topo_coords = np.zeros((3, self.nprocs), 
                                    order='F', dtype=np.int32)


    def setNumProcs(self):
        self.nprocs = self.npx * self.npy * self.npz

    def gridDims(self):
        return [self.npx, self.npy, self.npz]

    def createTopo(self, realm_comm, icoords=None):
        if self._topo_comm is None:
            self._topo_comm = realm_comm.Create_cart([self.npx, 
                                                     self.npy, self.npz])
            if icoords is None:
                self._createCartTopo()
            else:
                self._createCartTopoFromCoords()
            self.created = True
            self.checkProcs()

    def _createCartTopo(self):
        for rank in xrange(self.nprocs):
            cart_coords_l = self._topo_comm.Get_coords(rank)
            self._topo_coords[:, rank] = cart_coords_l 
        self._topo_coords += 1
    
    def _createCartTopoFromCoords(self):
        pass


class Packer(object):
    def __init__(self, name, region, callback, data_len):
        self.region = toCPLArray(region, np.int32)
        self.callback = callback
        self.name = name
        self.buffer = np.zeros((data_len))

    def __call__(self, context):
        self.callback(self.region, socket)



class TopoMap(object):
    def __init__(self, topology=None, grid=None):
        self.topo_map = None
        if (topology is not None) and (grid is not None):
            self.createMap(topology, grid)

    def createMap(self, topology, grid):
        self._create(topology, grid)
        nprocs = topology._topo_comm.Get_size()
        for rank in xrange(nprocs):
            self._computeProcBounds(rank, topology, grid)

    def _create(self, topology, grid):
        pass

    def _computeProcBounds(self, topology=None, grid=None):
        pass

class CartTopoMap(TopoMap):
    def __init__(self, topology=None, grid=None):
        TopoMap.__init__(self, topology, grid)

    def _computeProcBounds(self, rank, topology=None, grid=None):
        cart_coords = topology._topo_comm.Get_coords(rank)
        [x, y, z] = cart_coords
        self.topo_map[0, rank] = x*self.ncxl + 1
        self.topo_map[1, rank] = self.topo_map[0, rank] + self.ncxl - 1
        self.topo_map[2, rank] = y*self.ncyl + 1
        self.topo_map[3, rank] = self.topo_map[2, rank] + self.ncyl - 1
        self.topo_map[4, rank] = z*self.nczl + 1
        self.topo_map[5, rank] = self.topo_map[4, rank] + self.nczl - 1

    def _create(self, topology, grid):
        self.ncxl = grid.ncx / topology.npx
        self.ncyl = grid.ncy / topology.npy
        self.nczl = grid.ncz / topology.npz
        nprocs = topology._topo_comm.Get_size() 
        self.topo_map = np.zeros((6, nprocs), order='F', dtype=np.int32)

class CartRegion(object):
    def __init__(self, xlo, xhi, ylo, yhi, zlo, zhi):
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi
        self.ncxl = xhi - xlo
        self.ncyl = yhi - ylo
        self.nczl = zhi - zlo
    
    def isInside(self, reg):
        return  self.xlo <= reg.xlo and \
                self.ylo <= reg.ylo and \
                self.zlo <= reg.zlo and \
                self.xhi >= reg.xhi and \
                self.yhi >= reg.yhi and \
                self.zhi >= reg.zhi





class Grid(object):
    def __init__(self, xg=None, yg=None, zg=None):
        self.created = False
        if (xg is not None and yg is not None and zg is not None):
            self.xg = toCPLArray(xg, np.float64)
            self.yg = toCPLArray(yg, np.float64)
            self.zg = toCPLArray(yg, np.float64)
            self.created = True
        self.glob_bounds = []

    def createGrid():
        pass

    def computeGlobBounds(self):
        pass

class CartGrid(Grid):
    def __init__(self, ncxyz, lxyz, xg=None, yg=None, zg=None):
        Grid.__init__(self, xg, yg, zg)
        if (ncxyz) and (lxyz):
            self.ncx = ncxyz[0]
            self.ncy = ncxyz[1]
            self.ncz = ncxyz[2]

            if not self.created:
                print "Created grid"
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
            [self.xg, self.yg, self.zg] = create_CPL_cart_3Dgrid(self.ncx, 
                                                                self.ncy, 
                                                                self.ncz, 
                                                                dx, dy, dz)
            self.computeGlobBounds()
            self.created = True
        else:
            print "Exception not enough params"

    def computeGlobBounds(self):
        # Check if that is always true
        self.glob_bounds = np.array([[1, self.ncx], [1, self.ncy], [1, self.ncz]],
                                    order='F', dtype=np.int32)


            


class CartSocketMD(CPLSocket):
    _realm = CPL.MD_REALM
    def __init__(self, context=None):
        CPLSocket.__init__(self, context)
        self.initial_step = 0

    def _initialise(self, dt, npxyz, xyzL, topology): 
        CPLSocket._initialise(dt, npxyz, xyzL, topology)


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
        npxyz = toCPLArray(self.topology.gridDims(), np.int32)
        xyzL = toCPLArray([self.lx, self.ly, self.lz], np.float64)
        dummy_density = 1.0
        self.nsteps, self.initialstep = self._lib.md_init(self.dt, 
                                                self.topology._topo_comm, 
                                                self.topology._topo_coords, 
                                                npxyz, xyzL, dummy_density)

        rank = self.topology._topo_comm.Get_rank()
        my_coords = self.topology._topo_comm.Get_coords(rank)
        my_coords = toCPLArray(my_coords, dtype=np.int32)
        my_coords += 1# Fortran indices
        extents = self._lib.proc_extents(my_coords, CPL.MD_REALM)
        self.proc_region = CartRegion(extents[0], extents[1], extents[2],
                                    extents[3], extents[4], extents[5])

    def _initMD(self):
        pass

class CartSocketCFD(CPLSocket, CartTopoMap):
    _realm = CPL.CFD_REALM
    def __init__(self, context=None):
        CPLSocket.__init__(self, context)
        CartTopoMap.__init__(self)
        self.grid = None


    def _initialise(self, dt, npxyz, xyzL, topology, grid): 
        CPLSocket._initialise(dt, npxyz, xyzL, topology)
        self.grid = grid
        self.createMap(topology, grid)


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
        npxyz = toCPLArray([self.topology.npx, self.topology.npy, self.topology.npz], np.int32)
        xyzL = toCPLArray([self.lx, self.ly, self.lz], np.float64)
        ncxyz = toCPLArray([self.grid.ncx, self.grid.ncy, self.grid.ncz], np.int32)
        dummy_density = 1.0
        self.createMap(self.topology, self.grid)
        tm = self.topo_map
        # After slicing the array can be easily non-contiguous 
        iTmin = toCPLArray(tm[0, :])
        iTmax = toCPLArray(tm[1, :])
        jTmin = toCPLArray(tm[2, :])
        jTmax = toCPLArray(tm[3, :])
        kTmin = toCPLArray(tm[4, :])
        kTmax = toCPLArray(tm[5, :])
        if (self.topology._topo_comm.Get_rank() == 0):
            print "iTmin: " + str(iTmin)
            print "iTmax: " + str(iTmax)
            print "jTmin: " + str(jTmin)
            print "jTmax: " + str(jTmax)
            print "kTmin: " + str(kTmin)
            print "kTmax: " + str(kTmax)
            print "self.topo_map" + str(self.topo_map)
            print "npxyz: " + str(npxyz)
            print "xyzL: " + str(xyzL)
            print "ncxyz: " + str(ncxyz)
            print "topo_coords:" + str(self.topology._topo_coords)
            print "xg: " + str(self.grid.xg) + "len: " + str(len(self.grid.xg))
            print "yg: " + str(self.grid.yg)+ "len: " + str(len(self.grid.zg))
            print "zg: " + str(self.grid.zg)+ "len: " + str(len(self.grid.zg))
        ijkcmin = self.grid.glob_bounds[:, 0]
        ijkcmax = self.grid.glob_bounds[:, 1]
        self._lib.cfd_init(self.nsteps, self.dt, self.topology._topo_comm, 
                            self.topology._topo_coords, npxyz, xyzL, ncxyz,
                            dummy_density, ijkcmax, ijkcmin, iTmin, iTmax, 
                            jTmin, jTmax, kTmin, kTmax, self.grid.xg, 
                            self.grid.yg, self.grid.zg)

        rank = self.topology._topo_comm.Get_rank()
        self.proc_region = CartRegion(iTmin[rank], iTmax[rank], jTmin[rank],
                                    jTmax[rank], kTmin[rank], kTmax[rank])
                                        
    def _initCFD(self):
        pass






