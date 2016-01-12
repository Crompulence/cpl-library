from ctypes import c_int, c_double, c_bool, c_char_p, cdll, POINTER
import numpy as np
from numpy.ctypeslib import ndpointer
import os

__all__ = ["CPL"]

_libname = "cpl"
#TODO: Raise exception of library not loaded
_loaded = False


# All Python types except integers, strings, and unicode strings have to be 
# wrapped in their corresponding ctypes type, so that they can be converted 
# to the required C

_CPL_VARS = {"icmin_olap" : c_int,
            "jcmin_olap" : c_int,
            "kcmin_olap" : c_int}

class CPL:
    # Shared attribute containing the library
    _cpl_lib = cdll.LoadLibrary(os.path.abspath(_libname))

    def __init__(self):
        pass

    # py_test_python function 
    py_test_python = _cpl_lib.CPLC_test_python
    py_test_python.argtypes = \
            [c_int, c_double, c_bool, 
            ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'), 
            ndpointer(np.float64, ndim=2,  flags='aligned, f_contiguous'), 
            ndpointer(np.int32, shape=(2,),flags='aligned, f_contiguous'), 
            ndpointer(np.int32, shape=(2,),flags='aligned, f_contiguous'),]


    def test_python(self, int_p, doub_p, bool_p, int_pptr, doub_pptr):
        int_pptr_dims = np.array(int_pptr.shape, order='F',dtype=np.int32)
        doub_pptr_dims = np.array(doub_pptr.shape, order='F', dtype=np.int32)
        self.py_test_python(int_p, doub_p, bool_p, int_pptr, doub_pptr, 
                            int_pptr_dims, doub_pptr_dims)


    py_cfd_init = _cpl_lib.CPLC_cfd_init

    py_cfd_init.argtypes = \
            [c_int,                                                             #nsteps
            c_double,                                                           #dt
            c_int,                                                              #icomm_grid
            ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),         #icoord
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #npxyz_cfd
            ndpointer(np.float64, shape=(3), flags='aligned, f_contiguous'),    #xyzL
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ncxyz
            c_double,                                                           #density
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ijkcmax 
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ijkcmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #iTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #iTmax
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #jTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #jTmax
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #kTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #kTmax
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous'),       #xgrid
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous'),       #ygrid
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous')]       #zgrid


    def cfd_init(self, nsteps, dt, icomm_grid, icoord, npxyz_cfd, xyzL, ncxyz,
                density, ijkcmax, ijkcmin,iTmin, iTmax, jTmin, jTmax, kTmin, 
                kTmax, xgrid, ygrid, zgrid):

        self.py_cfd_init(nsteps, dt, icomm_grid, icoord, npxyz_cfd, xyzL, ncxyz,
                         density, ijkcmax, ijkcmin,iTmin, iTmax, jTmin, jTmax, 
                         kTmin, kTmax, xgrid, ygrid, zgrid)


    py_md_init = _cpl_lib.CPLC_md_init

    py_md_init.argtypes = \
            [POINTER(c_int),                                                             #nsteps
            POINTER(c_int),                                                              #icomm_grid
            ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),         #icoord
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #npxyz_cfd
            ndpointer(np.float64, shape=(3), flags='aligned, f_contiguous'),    #xyzL
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ncxyz
            c_double,                                                           #density
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ijkcmax 
            ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),     #ijkcmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #iTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #iTmax
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #jTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #jTmax
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #kTmin
            ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),         #kTmax
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous'),       #xgrid
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous'),       #ygrid
            ndpointer(np.float64, ndim=1, flags='aligned, f_contiguous')]       #zgrid


    def md_init(self, nstep, initialstep, dt, icom_grid, icoord, mpxyz_md, global_domain, density=1.0):
        self.py_md_init(nstep, initialstep, dt, icom_grid, icoord, mpxyz_md, global_domain, density)



    def get(self, var_name):
        try:
            var_type = _CPL_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_" + var_name)
        except KeyError:
            print "CPL-ERROR: CPL Library function '" + str(var_name) + "' not found!"
            #raise Function not found
        else:
            fun.restype = var_type
            return fun()

    py_gather = _cpl_lib.CPLC_gather
    py_gather.argtypes = \
                [ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'), 
                ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),
                ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'),
                ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'),
                ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous')]
    def gather(self, scatterarray, scatter_shape, limits, recvarray, recv_shape):
        pass

    py_scatter = _cpl_lib.CPLC_scatter
    py_scatter.argtypes = [ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'), ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'), ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous'), ndpointer(np.float64, ndim=2, flags='aligned, f_contiguous'), ndpointer(np.int32, ndim=2, flags='aligned, f_contiguous')]
    def scatter(self):
        pass



if __name__ == "__main__":
    _cpl_library = CPL(CPL_LIBRARY_PATH)
