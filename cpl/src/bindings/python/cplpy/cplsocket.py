from __future__ import division
from cpl import CPL, create_CPL_cart_3Dgrid, get_olap_limits, toCPLArray
import numpy as np

__all__ = ["CartSocketCFD", "CartSocketMD"]
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
        self.comm_method = CPL.SEND_RECEIVE
        # Flag to check if buffer is ready to be sent
        self.packed = False
        self.packers = {}
        self.unpackers = {}
        self.active_packer = None
        self.active_unpacker = None
        # Init domain and procs attributes
        self.lx = 0
        self.ly = 0
        self.lz = 0
        self.dt = 0.0
        self.topology = None
        self._createComms()
        self.proc_region = None
        self.olap_region = None

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
        self.packers[packer.name] = packer
        self.packers[packer.name].createBuffs(self)


    def registerUnpacker(self, unpacker):
        self.unpackers[unpacker.name] = unpacker
        self.unpackers[unpacker.name].createBuffs(self)


    
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
    def send(self, pack_name):
        self.active_packer = pack_name
        self.packers[self.active_packer].prepareBuffs(self)
        self.packers[self.active_packer](self)
        self._send()
    
    # Receive methods    
    def receive(self, unpack_name):
        self.active_unpacker = unpack_name
        self.unpackers[self.active_unpacker].prepareBuffs(self)
        self._receive()
        self.unpackers[self.active_unpacker](self)

    # Abstract methods
    def _init(self):
        pass

    def _send(self): 
        pass

    def _receive(self):
        pass

 
    def _computeProcRegion(self):
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



class CPLTopoMap(object):
    def __init__(self):
        self.topo_map = {}

    def createMap(self, topology, grid):
        self._create(topology, grid)
        nprocs = topology._topo_comm.Get_size()
        for rank in xrange(nprocs):
            self._computeProcBounds(rank, topology, grid)

    def _create(self, topology, grid):
        pass

    def _computeProcBounds(self, topology=None, grid=None):
        pass

class CartTopoMap(CPLTopoMap):
    def __init__(self, topology=None, grid=None):
        CPLTopoMap.__init__(self)

    def _computeProcBounds(self, rank, topology=None, grid=None):
        cart_coords = topology._topo_comm.Get_coords(rank)
        [x, y, z] = cart_coords
        self.topo_map["iTmin"][x] = x*self.ncxl + 1
        self.topo_map["iTmax"][x] = self.topo_map["iTmin"][x] + self.ncxl - 1
        self.topo_map["jTmin"][y] = y*self.ncyl + 1
        self.topo_map["jTmax"][y] = self.topo_map["jTmin"][y] + self.ncyl - 1
        self.topo_map["kTmin"][z] = z*self.nczl + 1
        self.topo_map["kTmax"][z] = self.topo_map["kTmin"][z] + self.nczl - 1

    def _create(self, topology, grid):
        self.ncxl = grid.ncx / topology.npx
        self.ncyl = grid.ncy / topology.npy
        self.nczl = grid.ncz / topology.npz
        self.topo_map["iTmin"] = np.zeros(topology.npx, order='F', dtype=np.int32)
        self.topo_map["iTmax"] = np.zeros(topology.npx, order='F', dtype=np.int32)
        self.topo_map["jTmin"] = np.zeros(topology.npy, order='F', dtype=np.int32)
        self.topo_map["jTmax"] = np.zeros(topology.npy, order='F', dtype=np.int32)
        self.topo_map["kTmin"] = np.zeros(topology.npz, order='F', dtype=np.int32)
        self.topo_map["kTmax"] = np.zeros(topology.npz, order='F', dtype=np.int32)

class CartRegion(object):
    def __init__(self, *args, **kwargs): 
        params_ok = True
        limits = kwargs.get("limits")
        if limits is not None:
            [self.xlo, self.xhi, self.ylo, 
            self.yhi, self.zlo, self.zhi] = kwargs["limits"]
        elif len(args) == 6:
            [self.xlo, self.xhi, self.ylo,
            self.yhi, self.zlo, self.zhi] = args
        else:
            params_ok = False
        if params_ok:
            self.ncxl = self.xhi - self.xlo + 1
            self.ncyl = self.yhi - self.ylo + 1
            self.nczl = self.zhi - self.zlo + 1
        else:
            print "Something is wrong with the parameters."
        
    def isInside(self, reg):
        return  self.xlo <= reg.xlo and \
                self.ylo <= reg.ylo and \
                self.zlo <= reg.zlo and \
                self.xhi >= reg.xhi and \
                self.yhi >= reg.yhi and \
                self.zhi >= reg.zhi

    def flatten(self):
        return toCPLArray([self.xlo, self.xhi, self.ylo, 
                        self.yhi, self.zlo, self.zhi], np.int32)

    def __repr__(self):
        return str(self.flatten())





class CPLGrid(object):
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

class CartGrid(CPLGrid):
    def __init__(self, ncxyz, lxyz, xg=None, yg=None, zg=None):
        CPLGrid.__init__(self, xg, yg, zg)
        if (ncxyz) and (lxyz):
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


    def _send(self):
        limits = self.packers[self.active_packer].region.flatten()
        if self.comm_method == CPL.GATHER_SCATTER:
            self._lib.gather(self.send_buff, limits, self.recv_buff)
        elif self.comm_method == CPL.SEND_RECEIVE:
            self._lib.send(self.send_buff, *limits)
        else:
            print "No recognised communication method"


    def _receive(self):
        limits = self.unpackers[self.active_unpacker].region.flatten()
        if self.comm_method == CPL.GATHER_SCATTER:
            self._lib.scatter(self.send_buff, limits, self.recv_buff)
        elif self.comm_method == CPL.SEND_RECEIVE:
            self._lib.recv(self.recv_buff, *limits)
        else:
            print "No recognised communication method"

    def _computeProcRegion(self):
        rank = self.topology._topo_comm.Get_rank()
        self.my_coords = self.topology._topo_comm.Get_coords(rank)
        self.my_coords = toCPLArray(self.my_coords, np.int32)
        self.my_coords += 1# Fortran indices
        extents = self._lib.proc_extents(self.my_coords, CPL.MD_REALM)
        self.proc_region = CartRegion(extents[0], extents[1], extents[2],
                                      extents[3], extents[4], extents[5])

    def _init(self):
        self._initMD()
        npxyz = toCPLArray(self.topology.gridDims(), np.int32)
        xyzL = toCPLArray([self.lx, self.ly, self.lz], np.float64)
        dummy_density = 1.0
        self.nsteps, self.initialstep = self._lib.md_init(self.dt, 
                                                self.topology._topo_comm, 
                                                self.topology._topo_coords, 
                                                npxyz, xyzL, dummy_density)
        self._computeProcRegion()
        self.olap_region = CartRegion(limits=get_olap_limits(self._lib))
    #    if (True):
    #        print str(self.topology._topo_comm.Get_rank())+" "+ "mycoords:" + str(self.my_coords) +" extents: " + str(self.proc_region.flatten())+" "+ str(self.topology._topo_comm)+ \
    #                        str(self.topology._topo_coords)+ str(npxyz)+ str(xyzL)
        

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


    def _send(self):
        limits = self.packers[self.active_packer].region.flatten()
        if self.comm_method == CPL.GATHER_SCATTER:
            self._lib.scatter(self.send_buff, limits, self.recv_buff)
        elif self.comm_method == CPL.SEND_RECEIVE:
            self._lib.send(self.send_buff, *limits)
        else:
            print "No recognised communication method"

    def _receive(self):
        limits = self.unpackers[self.active_unpacker].region.flatten()
        if self.comm_method == CPL.GATHER_SCATTER:
            self._lib.gather(self.send_buff, limits, self.recv_buff)
        elif self.comm_method == CPL.SEND_RECEIVE:
            self._lib.recv(self.recv_buff, *limits)
        else:
            print "No recognised communication method"


    def _computeProcRegion(self):
        rank = self.topology._topo_comm.Get_rank()
        my_coords = self.topology._topo_comm.Get_coords(rank)
        my_coords = toCPLArray(my_coords, np.int32)
        self.my_coords = my_coords
        my_coords += 1# Fortran indices
        extents = self._lib.proc_extents(my_coords, CPL.CFD_REALM)
        self.proc_region = CartRegion(extents[0], extents[1], extents[2],
                                    extents[3], extents[4], extents[5])


    def _initCFD(self):
        pass

    def _init(self):
        self._initCFD()
        self.createMap(self.topology, self.grid)
        npxyz = toCPLArray(self.topology.gridDims(), np.int32)
        xyzL = toCPLArray([self.lx, self.ly, self.lz], np.float64)
        ncxyz = toCPLArray([self.grid.ncx, self.grid.ncy, self.grid.ncz], np.int32)
        dummy_density = 1.0
        self.createMap(self.topology, self.grid)
        tm = self.topo_map
        # After slicing the array can easily be non-contiguous 
        iTmin = tm["iTmin"]
        iTmax = tm["iTmax"]
        jTmin = tm["jTmin"]
        jTmax = tm["jTmax"]
        kTmin = tm["kTmin"]
        kTmax = tm["kTmax"]
        ijkcmin = self.grid.glob_bounds[:, 0]
        ijkcmax = self.grid.glob_bounds[:, 1]
        self._lib.cfd_init(self.nsteps, self.dt, self.topology._topo_comm, 
                            self.topology._topo_coords, npxyz, xyzL, ncxyz,
                            dummy_density, ijkcmax, ijkcmin, iTmin, iTmax, 
                            jTmin, jTmax, kTmin, kTmax, self.grid.xg, 
                            self.grid.yg, self.grid.zg)

        self._computeProcRegion()
        self.olap_region = CartRegion(limits=get_olap_limits(self._lib))
#        if self.topology._topo_comm.Get_rank() == 2:
"""
        if (True):
            print str(self.topology._topo_comm.Get_rank())+" "+ "mycoords:" + str(self.my_coords+1) +" extents: " + str(self.proc_region.flatten())+" "+ str(self.nsteps)+ str(self.dt)+ str(self.topology._topo_comm)+ \
                            str(self.topology._topo_coords)+ str(npxyz)+ str(xyzL)+ str(ncxyz)+\
                            str(dummy_density)+ str(ijkcmax)+ str(ijkcmin)+ str(iTmin)+ str(iTmax)+\
                            str(jTmin)+ str(jTmax)+ str(kTmin)+ str(kTmax)
"""

