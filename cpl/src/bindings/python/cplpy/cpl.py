from ctypes import c_int, c_double, c_bool, c_void_p, byref, POINTER, util
import ctypes
from mpi4py import MPI
import numpy as np
from numpy.ctypeslib import ndpointer, load_library
import os


__all__ = ["CPL", "create_CPL_cart_3Dgrid", "cart_create",
           "get_olap_limits", "toCPLArray"]

class OpenMPI_Not_Supported(Exception):
    pass


# TODO: Raise exception of library not loaded
_loaded = False


# All Python types except integers, strings, and unicode strings have to be
# wrapped in their corresponding ctypes type, so that they can be converted
# to the required C

_CPL_GET_VARS = {"icmin_olap": c_int,
                 "jcmin_olap": c_int,
                 "kcmin_olap": c_int,
                 "icmax_olap": c_int,
                 "jcmax_olap": c_int,
                 "kcmax_olap": c_int,
                 "ncx": c_int,
                 "ncy": c_int,
                 "ncz": c_int,
                 "npx_md": c_int,
                 "npy_md": c_int,
                 "npz_md": c_int,
                 "npx_cfd": c_int,
                 "npy_cfd": c_int,
                 "npz_cfd": c_int,
                 "overlap": c_int,
                 "xl_md": c_double,
                 "yl_md": c_double,
                 "zl_md": c_double,
                 "xl_cfd": c_double,
                 "yl_cfd": c_double,
                 "zl_cfd": c_double}

_CPL_SET_VARS = {"output_mode": c_int}


class CPL:
    # Shared attribute containing the library
    CFD_REALM = 1
    MD_REALM = 2
    GATHER_SCATTER = 1
    SEND_RECEIVE = 2
    NULL_REALM = 0
    _libname = "libcpl"
    _lib_path = os.environ.get("CPL_LIBRARY_PATH")
    #Try using CPL_LIBRARY_PATH and if not look 
    #in system path (or ctype path)
    try:
        _cpl_lib = load_library(_libname, _lib_path)
    except OSError:
        print("WARNING - "+ _libname +" not found under path specified by CPL_LIBRARY_PATH variable: "+ _lib_path)
        print("Attempting to find system version")
        _lib_path = util.find_library(_libname)
    # If not found, try to look around the current directory system
    # for the library
    try:
        _cpl_lib = load_library(_libname, _lib_path)
    except OSError:
        print("WARNING - "+ _libname +" not found in system path or $CPL_LIBRARY_PATH: "+ _lib_path)
        print("Attempting to search current folder system")
        trydir = ''; Found=False
        for level in range(10):
            if not Found:
                trydir += "../"
            for root, dirs, files in os.walk(trydir):
                print(level,trydir,root)
                if "libcpl.so" in files:
                    print("A version of "+ _libname +" is found under path: "+ root + ". Attempting to use this..")
                    Found=True
                    _lib_path = root
                    break

        #Note this is a security issue, maybe better to inform the user here?
        if found:
            _cpl_lib = load_library(_libname, _lib_path)
        else:
            print(_libname +" not found, please ensure you have compiled correctly")
            print( " and set CPL_LIBRARY_PATH to its location" )

    def __init__(self):
        pass

    # py_test_python function
    py_test_python = _cpl_lib.CPLC_test_python
    py_test_python.argtypes = \
        [c_int,
         c_double,
         c_bool,
         ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),
         ndpointer(np.float64, ndim=2,  flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(2,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(2,), flags='aligned, f_contiguous')]

    def test_python(self, int_p, doub_p, bool_p, int_pptr, doub_pptr):
        int_pptr_dims = np.array(int_pptr.shape, order='F', dtype=np.int32)
        doub_pptr_dims = np.array(doub_pptr.shape, order='F', dtype=np.int32)
        self.py_test_python(int_p, doub_p, bool_p, int_pptr, doub_pptr,
                            int_pptr_dims, doub_pptr_dims)

    _py_create_comm = _cpl_lib.CPLC_create_comm
    #OpenMPI comm greater than c_int
    if MPI._sizeof(MPI.Comm) == ctypes.sizeof(c_int): 
        _py_create_comm.argtypes = [c_int, POINTER(c_int)]
    else:
        excptstr ="Problem is in create_comm wrapper, as the OpenMPI COMM handle is not "
        excptstr += "an integer, c_void_p should be used so C bindings needs something like **void" 
        excptstr += "(No idea what to do in the Fortran code, maybe MPI_COMM_f2C required)"
        raise OpenMPI_Not_Supported(excptstr)
        _py_create_comm.argtypes = [c_int, POINTER(c_void_p)]

    def create_comm(self, calling_realm):

        # build a communicator mpi4py python object from the
        # handle returned by the CPL_create_comm function.
        if MPI._sizeof(MPI.Comm) == ctypes.sizeof(c_int): 
            COMM_ctype = c_int 
        else: 
            COMM_ctype = c_void_p 

        #Call create comm
        returned_realm_comm = COMM_ctype()
        self._py_create_comm(calling_realm, byref(returned_realm_comm))

        # Use an intracomm object as the template and override value
        newcomm = MPI.Intracomm()
        newcomm_ptr = MPI._addressof(newcomm) 
        comm_val = COMM_ctype.from_address(newcomm_ptr)
        comm_val.value = returned_realm_comm.value

        return newcomm

    py_cfd_init = _cpl_lib.CPLC_cfd_init

    py_cfd_init.argtypes = \
        [c_int,
         c_double,
         c_int,
         ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_double,
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'),
         ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'),
         ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous')]

    def cfd_init(self, nsteps, dt, icomm_grid, icoord, npxyz_cfd, xyzL, ncxyz,
                 density, ijkcmax, ijkcmin, iTmin, iTmax, jTmin, jTmax, kTmin,
                 kTmax, xgrid, ygrid, zgrid):

        self.py_cfd_init(nsteps, dt, MPI._handleof(icomm_grid), icoord,
                         npxyz_cfd, xyzL, ncxyz, density, ijkcmax, ijkcmin,
                         iTmin, iTmax, jTmin, jTmax, kTmin, kTmax, xgrid,
                         ygrid, zgrid)

    py_md_init = _cpl_lib.CPLC_md_init

    py_md_init.argtypes = \
        [POINTER(c_int),
         POINTER(c_int),
         c_double,
         c_int,
         ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous'),
         c_double]

    def md_init(self, dt, icom_grid, icoord, npxyz_md, global_domain,
                density=1.0):
        nsteps = c_int()
        initialstep = c_int()
        self.py_md_init(byref(nsteps), byref(initialstep), dt,
                        MPI._handleof(icom_grid), icoord, npxyz_md,
                        global_domain, density)
        return (nsteps.value, initialstep.value)

    def get(self, var_name):
        try:
            var_type = _CPL_GET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_" + var_name)
        except KeyError:
            print "CPL-ERROR: CPL Library function '" + \
                  str(var_name) + "' not found!"
        else:
            fun.restype = var_type
            return fun()

    def set(self, var_name, value):
        try:
            var_type = _CPL_SET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_set_" + var_name)
        except KeyError:
            print "CPL-ERROR: CPL Library function '" + \
                  str(var_name) + "' not found!"
            # Raise Function not found
        else:
            fun.argtypes = [var_type]
            return fun(var_type(value))

    def _type_check(self, A):
        if type(A) is list:
            A = np.array(A)
        if not A.flags["F_CONTIGUOUS"]:
            A = np.require(A, requirements=['F'])
        if not A.flags["ALIGNED"]:
            A = np.require(A, requirements=['A'])
        #if (len(A.shape) != 4):
        #    raise 
        return A

    py_gather = _cpl_lib.CPLC_gather

    py_gather.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous')]

    def gather(self, gather_array, limits, recv_array):
        gather_array = self._type_check(gather_array)
        recv_array = self._type_check(recv_array)
        gather_shape = np.array(gather_array.shape, order='F', dtype=np.int32)
        recv_shape = np.array(recv_array.shape, order='F', dtype=np.int32)
        self.py_gather(gather_array, gather_shape, limits, recv_array,
                       recv_shape)

        return recv_array

    py_scatter = _cpl_lib.CPLC_scatter

    py_scatter.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous')]

    def scatter(self, scatter_array, limits, recv_array):
        scatter_array = self._type_check(scatter_array)
        recv_array = self._type_check(recv_array)
        scatter_shape = np.array(scatter_array.shape, order='F',
                                 dtype=np.int32)
        recv_shape = np.array(recv_array.shape, order='F', dtype=np.int32)
        self.py_scatter(scatter_array, scatter_shape, limits,
                        recv_array, recv_shape)
        return recv_array

    py_proc_extents = _cpl_lib.CPLC_proc_extents
    py_proc_extents.argtypes = \
        [ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_int,
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    def proc_extents(self, coord, realm):
        coord = self._type_check(coord)
        extents = np.zeros(6, order='F', dtype=np.int32)
        self.py_proc_extents(coord, realm, extents)
        return extents

    py_proc_portion = _cpl_lib.CPLC_proc_portion
    py_proc_portion.argtypes = \
        [ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_int,
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    def proc_portion(self, coord, realm, limits):
        coord = self._type_check(coord)
        portion = np.zeros(6, order='F', dtype=np.int32)
        self.py_proc_portion(coord, realm, limits, portion)
        return portion

    py_send = _cpl_lib.CPLC_send
    py_send.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         c_int,
         c_int, c_int,
         c_int, c_int,
         c_int, c_int, POINTER(c_bool)]

    def send(self, asend, icmin, icmax, jcmin, jcmax, kcmin, kcmax):
        asend = self._type_check(asend)
        asend_shape = np.array(asend.shape, order='F', dtype=np.int32)
        ndims = asend.ndim
        send_flag = c_bool()
        self.py_send(asend, asend_shape, ndims, icmin, icmax, jcmin, jcmax,
                     kcmin, kcmax, byref(send_flag))
        return send_flag.value

    py_recv = _cpl_lib.CPLC_recv
    py_recv.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         c_int,
         c_int, c_int,
         c_int, c_int,
         c_int, c_int, POINTER(c_bool)]

    def recv(self, arecv, icmin, icmax, jcmin, jcmax, kcmin, kcmax):
        arecv = self._type_check(arecv)
        arecv_shape = np.array(arecv.shape, order='F', dtype=np.int32)
        ndims = arecv.ndim
        recv_flag = c_bool()
        self.py_recv(arecv, arecv_shape, ndims, icmin, icmax, jcmin, jcmax,
                     kcmin, kcmax, byref(recv_flag))
        return arecv, recv_flag.value




def create_CPL_cart_3Dgrid(ncx, ncy, ncz, dx, dy, dz):
    xg = np.zeros((ncx + 1, ncy + 1), order='F', dtype=np.float64)
    yg = np.zeros((ncx + 1, ncy + 1), order='F', dtype=np.float64)
    zg = np.zeros(ncz + 1, order='F', dtype=np.float64)
    for i in xrange(ncx + 1):
        for j in xrange(ncy + 1):
            xg[i, j] = i*dx
            yg[i, j] = j*dy
    for k in xrange(ncz + 1):
        zg[k] = k*dz
    return (xg, yg, zg)


def cart_create(old_comm, dims, periods, coords):
    dummy_cart_comm = old_comm.Create_cart(dims, periods)
    temp_comm = old_comm.Split(0, dummy_cart_comm.Get_rank())
    new_cart_comm = temp_comm.Create_cart(dims, periods)
    comm_coords = new_cart_comm.Get_coords(new_cart_comm.Get_rank())
    if (not (coords == comm_coords).all()):
            print "Not good!"
            exit()
    return new_cart_comm


def get_olap_limits(lib):
    olap_limits = np.zeros(6, order='F', dtype=np.int32)
    olap_limits[0] = lib.get("icmin_olap")
    olap_limits[1] = lib.get("icmax_olap")
    olap_limits[2] = lib.get("jcmin_olap")
    olap_limits[3] = lib.get("jcmax_olap")
    olap_limits[4] = lib.get("kcmin_olap")
    olap_limits[5] = lib.get("kcmax_olap")
    return olap_limits


def toCPLArray(arr, arr_type=None):
    if type(arr) == np.ndarray:
        if not arr.flags["F_CONTIGUOUS"]:
            return np.asfortranarray(arr, dtype=arr.dtype)
    else:
        if arr_type is not None:
            return np.asfortranarray(arr, dtype=arr_type)
        else:
            print "Non-numpy arrays require argument arr_type"

if __name__ == "__main__":
    lib = CPL()
