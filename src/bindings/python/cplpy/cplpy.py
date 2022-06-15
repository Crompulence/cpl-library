
from ctypes import c_char_p, c_char, c_int, c_double, c_bool, c_void_p, byref, POINTER, util, pointer, cdll
import ctypes
import mpi4py
from distutils.version import StrictVersion
from mpi4py import MPI
import numpy as np
from numpy.ctypeslib import ndpointer, load_library
from functools import wraps
import os
import time
import subprocess as sp
from subprocess import STDOUT, check_output, CalledProcessError
import shutil
import pickle
import errno


__all__ = ["CPL", "cart_create", "run_test", "prepare_config", "parametrize_file"]


#class OpenMPI_Not_Supported(Exception):
#    pass

class CPLLibraryNotFound(Exception):
    pass

class mpi4py_version_error(Exception):
    pass


# TODO: Raise exception of library not loaded
_loaded = False

# All Python types except integers, strings, and unicode strings have to be
# wrapped in their corresponding ctypes type, so that they can be converted
# to the required C

_CPL_GET_VARS = {"icmin_olap": c_int, "jcmin_olap": c_int, "kcmin_olap": c_int,
                 "icmax_olap": c_int, "jcmax_olap": c_int, "kcmax_olap": c_int,
                 "icmin_cnst": c_int, "jcmin_cnst": c_int, "kcmin_cnst": c_int,
                 "icmax_cnst": c_int, "jcmax_cnst": c_int, "kcmax_cnst": c_int,
                 "icmin_bnry": c_int, "jcmin_bnry": c_int, "kcmin_bnry": c_int,
                 "icmax_bnry": c_int, "jcmax_bnry": c_int, "kcmax_bnry": c_int,
                 "ncx": c_int, "ncy": c_int, "ncz": c_int,
                 "npx_md": c_int, "npy_md": c_int, "npz_md": c_int,
                 "npx_cfd": c_int, "npy_cfd": c_int, "npz_cfd": c_int,
                 "overlap": c_int, "xl_md": c_double, "yl_md": c_double,
                 "nsteps_md": c_int, "nsteps_cfd": c_int, "nsteps_coupled": c_int,
                 "zl_md": c_double, "xl_cfd": c_double, "yl_cfd": c_double,
                 "zl_cfd": c_double, "dx" : c_double, "dy" : c_double, "dz" : c_double,
                 "x_orig_cfd": c_double,"y_orig_cfd": c_double,"z_orig_cfd": c_double,
                 "x_orig_md": c_double,"y_orig_md": c_double,"z_orig_md": c_double,
                 "timestep_ratio": c_int
                 }

_CPL_SET_VARS = {"output_mode": c_int}

_CPL_FILE_VARS_TYPES = {"CPL_INT", "CPL_DOUBLE", "CPL_INT_ARRAY", "CPL_DOUBLE_ARRAY"}

class CPL_VAR_TYPES():
    INT = 1
    DOUBLE = 2
    BOOL = 3
    STRING = 4
    INT_ARRAY = 5
    DOUBLE_ARRAY = 6
    BOOL_ARRAY = 7
    STRING_ARRAY = 8


_CPL_GET_FILE_VARS = {CPL_VAR_TYPES.DOUBLE: ("get_real_param", c_double), 
                      CPL_VAR_TYPES.DOUBLE_ARRAY:("get_real_array_param", POINTER(c_double)),
                      CPL_VAR_TYPES.INT: ("get_int_param", c_int), 
                      CPL_VAR_TYPES.INT_ARRAY: ("get_int_array_param", POINTER(c_int)),
                      CPL_VAR_TYPES.BOOL: ("get_boolean_param", c_bool), 
                      CPL_VAR_TYPES.BOOL_ARRAY: ("get_boolean_array_param", POINTER(c_bool)),
                      CPL_VAR_TYPES.STRING: ("get_string_param", c_char_p), 
                      CPL_VAR_TYPES.STRING_ARRAY: ("get_string_array_param", POINTER(c_char_p))}
# Decorator to abort all processes if an exception is thrown. This
# avoids getting blocked when the exception do not occurs in every
# process.
def abortMPI(func):
    @wraps(func)
    def handleExcepts(self, *args, **kwargs):
        retval = None
        try:
            retval = func(self, *args, **kwargs)
        except Exception as e:
            print (e)
            # Dirty workaround to let the output be printed.
            time.sleep(2)
            MPI.COMM_WORLD.Abort(errorcode=1)
        else:
            return retval
    return handleExcepts


class CPL:
    # Shared attribute containing the library
    CFD_REALM = 1
    MD_REALM = 2
    GATHER_SCATTER = 1
    SEND_RECEIVE = 2
    NULL_REALM = 0
    _libname = "libcpl"
    try:
        _lib_path = os.environ["CPL_LIBRARY_PATH"]
        if os.path.exists(_lib_path + "/"+ _libname + ".so"):
            _cpl_lib = load_library(_libname, _lib_path)
        else:
            raise CPLLibraryNotFound("Compiled CPL library libcpl.so not found at " + _lib_path + "/"+ _libname + ".so")
    except KeyError as e:
        print(("CPL info: ", "CPL_LIBRARY_PATH not defined. Looking in system directories..."))
        try:
            _cpl_lib = cdll.LoadLibrary(_libname + ".so")
            print(("CPL info: ", "Success!"))
        except OSError as e:
            raise CPLLibraryNotFound("Library libcpl.so not found!")
            #TODO: Check this
            #time.sleep(2)
            #MPI.COMM_WORLD.Abort(errorcode=1)

    # Check for JSON support by cheking if load_param_file symbol exists
    JSON_SUPPORT = True
    try:
        _cpl_lib.CPLC_load_param_file
    except:
        JSON_SUPPORT = False

    def __init__(self):
        self._var = POINTER(POINTER(c_char_p))
        self.realm = None

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

    @abortMPI
    def test_python(self, int_p, doub_p, bool_p, int_pptr, doub_pptr):
        int_pptr_dims = np.array(int_pptr.shape, order='F', dtype=np.int32)
        doub_pptr_dims = np.array(doub_pptr.shape, order='F', dtype=np.int32)
        self.py_test_python(int_p, doub_p, bool_p, int_pptr, doub_pptr,
                            int_pptr_dims, doub_pptr_dims)

    #NOTE: Using CPLC_init_Fort and Comm.f2py() and Comm.py2f() we achieve integration
    #      with MPICH and OpenMPI seamlessly. mpi4py >= 2.0.0 is needed.
    if StrictVersion(mpi4py.__version__) < StrictVersion('2.0.0'):
        raise mpi4py_version_error("Comm.f2py() and Comm.py2f()" + 
                                   " require mpi4py >= 2.0.0")

    #Detect if OpenMPI or MPICH
    mpicc=str(mpi4py.get_config()['mpicc'])
    mpishow=sp.check_output(["mpicc","-show"]).decode("utf-8")
    if ("open" in mpicc or "open" in mpishow):
        MPI_version = "OPENMPI"
        ompi_info = sp.check_output("ompi_info").split("\n")
        for m in ompi_info:
            if ("Open MPI:" in m):
                ompi_major_version_no = int(m.split(":")[-1].split(".")[0])            
    elif ("mpich" in mpicc or "mpich" in mpishow):
        MPI_version = "MPICH"      
    else:
        print(("UNKNOWN MPI VERSION FROM ", mpicc))
        MPI_version = "UNKNOWN"

    _py_init = _cpl_lib.CPLC_init_Fort
    _py_init.argtypes = [c_int, POINTER(c_int)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def init(self, calling_realm):
        """
            (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively 
            and create intercommunicator between CFD and MD
            
            **Remarks**
            
            Assumes MPI has been initialised `MPI_init` and communicator MPI_COMM_WORLD exists
            and contains all processors in both CFD and MD regions
            
            
            **Synopsis**
            
            .. code-block:: python
            
              CPL.init(callingrealm)    
            
            **Inputs**
            
             - *callingrealm*
            
               - Should identify calling processor as either CFD_REALM (integer with value 1) or MD_REALM (integer with value 2).
             
            
            **Outputs**
            
             - RETURNED_REALM_COMM 
            
               - Communicator based on callingrealm value local to CFD or MD processor and resulting from the split of MPI_COMM_WORLD
            
            **Example**
            
            .. literalinclude:: ../../../examples/cpl_init/md_init.py           
            
            **Errors**
            
                COUPLER_ERROR_REALM  = 1                    wrong realm value
                COUPLER_ERROR_ONE_REALM = 2                 one realm missing
                COUPLER_ERROR_INIT = 3         ! initialisation error
        """
        self.realm = calling_realm
        MPI_Comm_int = c_int()
        self._py_init(calling_realm, byref(MPI_Comm_int))
        self.COMM = MPI.Comm.f2py(MPI_Comm_int.value)
        return self.COMM 


    if JSON_SUPPORT:
        _py_load_param_file = _cpl_lib.CPLC_load_param_file
        _py_load_param_file.argtypes = [c_char_p]

    @abortMPI
    def load_param_file(self, fname):
        self._py_load_param_file(c_char_p(fname), c_int(len(fname)))


    if JSON_SUPPORT:
        _py_close_param_file = _cpl_lib.CPLC_close_param_file

    @abortMPI
    def close_param_file(self):
        self._py_close_param_file()


    @abortMPI
    def get_file_var(self, section, var_name, var_type):
        try:
            fun_name = _CPL_GET_FILE_VARS[var_type][0]
            var_ctype = _CPL_GET_FILE_VARS[var_type][1]

            fun = getattr(self._cpl_lib, "CPLC_" + fun_name)
            fun.argtypes =  [c_char_p,
                             c_char_p,
                             POINTER(var_ctype)]

        except KeyError:
            print(("CPL-ERROR: CPL Library function '" +
                   str(fun_name) + "' not found!"))
            raise KeyError
        else:
            self._var = var_ctype()

            if ("array" in fun_name):
                print ("ENTRO")
                var_len = c_int()
                fun.argtypes.append(POINTER(c_int))
                print ("EY")
                fun(c_char_p(section), c_char_p(var_name), byref(self._var), byref(var_len))
                print(("len:" , var_len.value))
                #print (self._var[0])
                #print (byref(var[0]))
                a = ([self._var[i] for i in range(var_len.value)])
                return a
            else:
                fun(c_char_p(section), c_char_p(var_name), byref(self._var))
                return self._var.value
    _py_finalize = _cpl_lib.CPLC_finalize

    @abortMPI
    def finalize(self):
        self._py_finalize()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_setup_cfd = _cpl_lib.CPLC_setup_cfd_Fort

    py_setup_cfd.argtypes = \
        [c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def setup_cfd(self, icomm_grid, xyzL, 
                        xyz_orig, ncxyz):
        """
            Initialisation routine for coupler module - Every variable is sent and stored
            to ensure both md and cfd region have an identical list of parameters
           
            **Remarks**
            
            Assumes CPL has been initialised `CPL.init` and communicator MD_REALM exists

            **Synopsis**
            
            .. code-block:: python
            
              CPL.setup_cfd(icomm_grid, xyzL, xyz_orig, ncxyz)
            
            **Inputs**
            
             - *icomm_grid*
            
               - Communicator based on CFD processor topology returned from a call to MPI_CART_CREATE.
             - *xyzL*
            
               - CFD domain size.
             - *xyz_orig*
            
               - CFD origin.
             - *ncxyz*
            
               - Number of CFD cells in global domain.

        """

        if (  ((type(icomm_grid) is list) and (len(icomm_grid) is 3))
           or ((type(icomm_grid) is np.array) and (icomm_grid.shape[0] is 3))):
            icomm_grid = self.COMM.Create_cart([icomm_grid[0], 
                                                icomm_grid[1], 
                                                icomm_grid[2]])

        if ((type(xyzL) is list) or 
            (xyzL.dtype != np.float64) or 
            (not xyzL.flags["F_CONTIGUOUS"])):
            xyzL = np.array(xyzL, order='F', dtype=np.float64)

        if ((type(xyz_orig) is list) or 
            (xyz_orig.dtype != np.float64) or 
            (not xyz_orig.flags["F_CONTIGUOUS"])):
            xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)

        if ((type(ncxyz) is list) or 
            (ncxyz.dtype != np.int32) or 
            (not ncxyz.flags["F_CONTIGUOUS"])):
            ncxyz = np.array(ncxyz, order='F', dtype=np.int32)

        self.py_setup_cfd(icomm_grid.py2f(), xyzL,
                          xyz_orig, ncxyz)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


    py_setup_md = _cpl_lib.CPLC_setup_md_Fort

    py_setup_md.argtypes = \
        [c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def setup_md(self, icomm_grid, xyzL, xyz_orig):
        """
            Initialisation routine for coupler module - Every variable is sent and stored
            to ensure both md and cfd region have an identical list of parameters
            
            **Remarks**
            
            Assumes CPL has been initialised `CPL.init` and communicator MD_REALM exists
            
            **Synopsis**
            
            .. code-block:: python
            
              CPL.md_init(icomm_grid, xyzL, xyz_orig)
            
            **Inputs**
            
             - *icomm_grid*
            
               - Communicator based on MD processor topology returned from a call to MPI_CART_CREATE.
             - *xyzL*
            
               - MD domain size.
             - *xyz_orig*
            
               - MD origin.

        """
        if (  ((type(icomm_grid) is list) and (len(icomm_grid) is 3))
           or ((type(icomm_grid) is np.array) and (icomm_grid.shape[0] is 3))):
            icomm_grid = self.COMM.Create_cart([icomm_grid[0], 
                                                icomm_grid[1], 
                                                icomm_grid[2]])

        if ((type(xyzL) is list) or 
            (xyzL.dtype != np.float64) or 
            (not xyzL.flags["F_CONTIGUOUS"])):
            xyzL = np.array(xyzL, order='F', dtype=np.float64)

        if ((type(xyz_orig) is list) or 
            (xyz_orig.dtype != np.float64) or 
            (not xyz_orig.flags["F_CONTIGUOUS"])):
            xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)

        self.py_setup_md(icomm_grid.py2f(), xyzL, xyz_orig)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_proc_extents = _cpl_lib.CPLC_proc_extents
    py_proc_extents.argtypes = \
        [ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_int,
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def proc_extents(self, coord, realm):
        coord = self._type_check(coord)
        extents = np.zeros(6, order='F', dtype=np.int32)
        self.py_proc_extents(coord, realm, extents)
        return extents

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_my_proc_extents = _cpl_lib.CPLC_my_proc_extents
    py_my_proc_extents.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def my_proc_extents(self):
        extents = np.zeros(6, order='F', dtype=np.int32)
        self.py_my_proc_extents(extents)
        return extents

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_proc_portion = _cpl_lib.CPLC_proc_portion
    py_proc_portion.argtypes = \
        [ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_int,
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def proc_portion(self, coord, realm, limits):
        """
            Get maximum and minimum cell indices, i.e. the 'portion', of the
            input cell extents 'limits' that is contributed by the processor
            specified by processor coord. 
        
            **Remarks**
        
            Assumes the coupler has been initialised with `CPL.init <#cplpy.CPL.init>`_ and 
            topological mapping has been setup using either `CPL.setup_md <#cplpy.CPL.setup_md>`_ 
            or `CPL.setup_cfd <#cplpy.CPL.setup_cfd>`_ as appropriate.
            - Note: limits(6) and portion(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
        
            **Synopsis**
            
            .. code-block:: python
        
              cpl.proc_portion(coord, realm, limits)
        
            **Inputs**
            
            - coord
            
                - processor cartesian coordinate, list or numpy array of 3 integers 
            - realm
            
                - cfd_realm (1) or md_realm (2) (integer) 
            - limits
            
                - Array of cell extents that specify the input region, list or numpy array of 6 integers  
        
        
            **Outputs**
        
            - portion
                - Array of cell extents that define the local processor's
                  contribution to the input region 'limits', numpy array of 6 integers  
        
        """

        coord = self._type_check(coord)
        limits = self._type_check(limits)
        portion = np.zeros(6, order='F', dtype=np.int32)
        self.py_proc_portion(coord, realm, limits, portion)
        return portion

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_my_proc_portion = _cpl_lib.CPLC_my_proc_portion
    py_my_proc_portion.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def my_proc_portion(self, limits):
        """
            Get maximum and minimum cell indices, i.e. the 'portion' on calling process. 
        
            **Remarks**
        
            Assumes the coupler has been initialised with `CPL.init <#cplpy.CPL.init>`_ and 
            topological mapping has been setup using either `CPL.setup_md <#cplpy.CPL.setup_md>`_ 
            or `CPL.setup_cfd <#cplpy.CPL.setup_cfd>`_ as appropriate.
            - Note: limits(6) and portion(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
        
            **Synopsis**
            
            .. code-block:: python
        
                CPL.my_proc_portion(limits)
        
            **Inputs**
            
            - limits
            
                - Array of cell extents that specify the input region, list or numpy array of 6 integers  
        
        
            **Outputs**
        
            - portion
                - Array of cell extents that define the local processor's
                  contribution to the input region 'limits', numpy array of 6 integers  
        
        """
        limits = self._type_check(limits)
        portion = np.zeros(6, order='F', dtype=np.int32)
        self.py_my_proc_portion(limits, portion)
        return portion

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_map_cfd2md_coord = _cpl_lib.CPLC_map_cfd2md_coord
    py_map_cfd2md_coord.argtypes = \
        [ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_cfd2md_coord(self, coord_cfd):
        coord_cfd = self._type_check(coord_cfd)
        coord_md = np.zeros(3, order='F', dtype=np.float64)
        self.py_map_cfd2md_coord(coord_cfd, coord_md)
        return coord_md

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_map_md2cfd_coord = _cpl_lib.CPLC_map_md2cfd_coord
    py_map_md2cfd_coord.argtypes = \
        [ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_md2cfd_coord(self, coord_md):
        coord_md = self._type_check(coord_md)
        coord_cfd = np.zeros(3, order='F', dtype=np.float64)
        self.py_map_md2cfd_coord(coord_md, coord_cfd)
        return coord_cfd

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_map_glob2loc_cell = _cpl_lib.CPLC_map_glob2loc_cell
    py_map_glob2loc_cell.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_glob2loc_cell(self, limits, glob_cell):
        limits = self._type_check(limits)
        glob_cell = self._type_check(glob_cell)
        loc_cell = np.zeros(3, order='F', dtype=np.int32)
        self.py_map_glob2loc_cell(limits, glob_cell, loc_cell)
        return loc_cell

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_map_cell2coord = _cpl_lib.CPLC_map_cell2coord
    py_map_cell2coord.argtypes = \
        [c_int, c_int, c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_cell2coord(self, i, j, k):
        coord = np.zeros(3, order='F', dtype=np.float64)
        self.py_map_cell2coord(i, j, k, coord)
        return coord

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_map_coord2cell = _cpl_lib.CPLC_map_coord2cell
    py_map_coord2cell.argtypes = \
        [c_double, c_double, c_double,
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_coord2cell(self, x, y, z):
        cell = np.zeros(3, order='F', dtype=np.int32)
        self.py_map_coord2cell(x, y, z, cell)
        return cell

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_get_no_cells = _cpl_lib.CPLC_get_no_cells
    py_get_no_cells.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_no_cells(self, limits):
        limits = self._type_check(limits)
        no_cells = np.zeros(3, order='F', dtype=np.int32)
        self.py_get_no_cells(limits, no_cells)
        return no_cells

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    #Limits of overlap region
    py_get_olap_limits = _cpl_lib.CPLC_get_olap_limits
    py_get_olap_limits.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_olap_limits(self):
        limits = np.zeros(6, order='F', dtype=np.int32)
        self.py_get_olap_limits(limits)
        return limits

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    #Limits of contraint region
    py_get_cnst_limits = _cpl_lib.CPLC_get_cnst_limits
    py_get_cnst_limits.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_cnst_limits(self):
        limits = np.zeros(6, order='F', dtype=np.int32)
        self.py_get_cnst_limits(limits)
        return limits

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    #Limits of boundary region
    py_get_bnry_limits = _cpl_lib.CPLC_get_bnry_limits
    py_get_bnry_limits.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_bnry_limits(self):
        limits = np.zeros(6, order='F', dtype=np.int32)
        self.py_get_bnry_limits(limits)
        return limits

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_set_timing = _cpl_lib.CPLC_set_timing
    py_set_timing.argtypes = \
        [c_int, c_int, c_double]

    #Don't call abortMPI so it can be handled nicely in Python.
    #@abortMPI
    def set_timing(self, initialstep, nsteps, dt):
        class DepricatedException(Exception):
            """Raise Error as function should not be used"""
        raise DepricatedException("CPL set_timing is depricated and should not be used")
        self.py_set_timing(initialstep, nsteps, dt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_send = _cpl_lib.CPLC_send
    py_send.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'), 
         POINTER(c_bool)]

    py_send_min = _cpl_lib.CPLC_send_min
    py_send_min.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         POINTER(c_bool)]

    @abortMPI
    def send(self, asend, limits=None):
        """
            
            Send four dimensional array *asend* of data from all processors in the 
            current realm with data between global cell array *limits* to the 
            corresponding processors from the other realm.
            
            **Remarks**
            
            Assumes the coupler has been initialised with `CPL.init <#cplpy.CPL.init>`_ and 
            topological mapping has been setup using either `CPL.setup_md <#cplpy.CPL.setup_md>`_ 
            or `CPL.setup_cfd <#cplpy.CPL.setup_cfd>`_ as appropriate.
            
            **Synopsis**
            
            .. code-block:: python
            
              CPL.send(asend, limits=None)    
            
            **Inputs**
            
             - asend
            
               - Array of data to send. Should be a four dimensional Numpy array allocated using the number of cells on the current processor between the limits. For example, if overlap limits are 8 cells, between cells 0 and 7 split over 2 procesors, the first processor will have from 0 to 3 and the second from 4 to 7. This should be be obtained from `CPL.my_proc_portion(limits, portion) <#cplpy.CPL.my_proc_portion>`_ to allocate a Numpy array, or allocated using the helper function `CPL.get_arrays <#cplpy.CPL.get_arrays>`_ with appropriate sizes.
                .
             - limits [Optional]
            
               - Optional arguments limits specify if global limits of overlap region not used. These are in the global cell coordinates, and must match the corresponding recieve.
            
            **Outputs**
            
             - send_flag
            
               - Returned flag which indicates success or failure of send process
            
            **Example**

            This example links with the `CPL.recv <#cplpy.CPL.recv>`_ examples 
                        
            .. literalinclude:: ../../../examples/sendrecv_globcell/python/cfd_send_cells.py


        """
        asend = self._type_check(asend)
        asend_shape = np.array(asend.shape, order='F', dtype=np.int32)
        send_flag = c_bool()

        if limits is None:
            self.py_send_min(asend, asend_shape, byref(send_flag))
        else:
            self.py_send(asend, asend_shape, limits, byref(send_flag))

        return send_flag.value

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_recv = _cpl_lib.CPLC_recv
    py_recv.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'), 
         POINTER(c_bool)]

    py_recv_min = _cpl_lib.CPLC_recv_min
    py_recv_min.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         POINTER(c_bool)]


    @abortMPI
    def recv(self, arecv, limits=None):
        """
            
            Receive data from to local grid from the associated ranks from the other realm
            
            **Remarks**
            
            Assumes the coupler has been initialised with `CPL.init <#cplpy.CPL.init>`_ and 
            topological mapping has been setup using either `CPL.setup_md <#cplpy.CPL.setup_md>`_ 
            or `CPL.setup_cfd <#cplpy.CPL.setup_cfd>`_ as appropriate.
            
            **Synopsis**
            
            .. code-block:: python
            
              CPL.recv(arecv, limits=None)    

            **Inputs**
            
             - arecv
            
               - Array of data to recv. Should be a four dimensional Numpy array allocated using the number of cells on the current processor between the limits. For example, if overlap limits are 8 cells, between cells 0 and 7 split over 2 procesors, the first processor will have from 0 to 3 and the second from 4 to 7. This should be be obtained from `CPL.my_proc_portion(limits, portion) <#cplpy.CPL.my_proc_portion>`_ to allocate a Numpy array, or allocated using the helper function `CPL.get_arrays <#cplpy.CPL.get_arrays>`_ with appropriate sizes.

             - limits [Optional]
            
               - Limits in global cell coordinates, must be the same as corresponding send command.
            
            **Outputs**
            
             - recv_flag
            
               - Returned flag which indicates success or failure of recv process
            
            **Example**

            This example links with the `CPL.send <#cplpy.CPL.send>`_ examples 

            .. literalinclude:: ../../../examples/sendrecv_globcell/python/md_recv_cells.py

        """

        arecv = self._type_check(arecv)
        arecv_shape = np.array(arecv.shape, order='F', dtype=np.int32)
        recv_flag = c_bool()
        if limits is None:
            self.py_recv_min(arecv, arecv_shape, byref(recv_flag))
        else:
            self.py_recv(arecv, arecv_shape, limits, byref(recv_flag))
        return arecv, recv_flag.value

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_gather = _cpl_lib.CPLC_gather
    py_gather.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous')]

    @abortMPI
    def gather(self, gather_array, limits, recv_array):
        gather_array = self._type_check(gather_array)
        recv_array = self._type_check(recv_array)
        gather_shape = np.array(gather_array.shape, order='F', dtype=np.int32)
        recv_shape = np.array(recv_array.shape, order='F', dtype=np.int32)
        self.py_gather(gather_array, gather_shape, limits, recv_array,
                       recv_shape)

        return recv_array

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_scatter = _cpl_lib.CPLC_scatter
    py_scatter.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous')]

    @abortMPI
    def scatter(self, scatter_array, limits, recv_array):
        scatter_array = self._type_check(scatter_array)
        recv_array = self._type_check(recv_array)
        scatter_shape = np.array(scatter_array.shape, order='F',
                                 dtype=np.int32)
        recv_shape = np.array(recv_array.shape, order='F', dtype=np.int32)
        self.py_scatter(scatter_array, scatter_shape, limits,
                        recv_array, recv_shape)
        return recv_array

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_swaphalos = _cpl_lib.CPLC_swaphalos
    py_swaphalos.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous')]

    @abortMPI
    def swaphalos(self, A):
        A = self._type_check(A)
        A_shape = np.array(A.shape, order='F', dtype=np.int32)
        self.py_swaphalos(A, A_shape)
        return A

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_overlap = _cpl_lib.CPLC_overlap
    py_overlap.argtypes = []

    @abortMPI
    def overlap(self):
        self.py_overlap.restype = c_bool
        return self.py_overlap()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    py_is_proc_inside = _cpl_lib.CPLC_is_proc_inside
    py_is_proc_inside.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def is_proc_inside(self, region):
        self.py_is_proc_inside.restype = c_bool
        region = self._type_check(region)
        return self.py_is_proc_inside(region)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def get(self, var_name):
        try:
            var_type = _CPL_GET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_" + var_name)
        except KeyError:
            print(("CPL-ERROR: CPL Library function '" +
                   str(var_name) + "' not found!"))
            print ("Available options include: ")
            for var in _CPL_GET_VARS:
                print(var)
            raise KeyError
        else:
            fun.restype = var_type
            return fun()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def set(self, var_name, value):
        try:
            var_type = _CPL_SET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_set_" + var_name)
        except KeyError:
            print(("CPL-ERROR: CPL Library function '" +
                   str(var_name) + "' not found!"))
            raise KeyError
        else:
            fun.argtypes = [var_type]
            return fun(var_type(value))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def _type_check(self, A):
        if type(A) is list:
            ndtype = type(A[0])
            if ndtype == float:
                ndtype = np.float64
            elif ndtype == int:
                ndtype = np.int32
            A = np.asfortranarray(A, dtype=ndtype)
        if not A.flags["F_CONTIGUOUS"]:
            A = np.require(A, requirements=['F'])
        if not A.flags["ALIGNED"]:
            A = np.require(A, requirements=['A'])
        return A

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def get_arrays(self, recv_size, send_size):

        """
          Return recv array and send array based
          on constraint/boundary sizes

          **Example**

            A minimal example is possible using CPL.get_arrays, which shows paired sending and recv commands

            .. literalinclude:: ../../../examples/minimal_send_recv_mocks/minimal_MD.py
            .. literalinclude:: ../../../examples/minimal_send_recv_mocks/minimal_CFD.py

        """
        #Get constraint region
        cnst_limits = self.get_cnst_limits();
        cnst_portion = self.my_proc_portion(cnst_limits)
        cnst_ncxl, cnst_ncyl, cnst_nczl = self.get_no_cells(cnst_portion)

        #Get overlap region
        BC_limits = self.get_bnry_limits()
        BC_portion = self.my_proc_portion(BC_limits)
        BC_ncxl, BC_ncyl, BC_nczl = self.get_no_cells(BC_portion)

        #Allocate send and recv arrays
        recv_array = np.zeros((recv_size, BC_ncxl, BC_ncyl, BC_nczl), order='F', dtype=np.float64)
        send_array = np.zeros((send_size, cnst_ncxl, cnst_ncyl, cnst_nczl), order='F', dtype=np.float64)

        return recv_array, send_array

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    @abortMPI
    def dump_region(self, region, array, fname, comm, components={}, coords="mine"):
            lines = ""
            portion = self.my_proc_portion(region)
            cell_coords = np.array(3)
            dx = self.get("dx")
            dy = self.get("dy")
            dz = self.get("dz")
            def_func = lambda x : x
            components_dic = components
            if callable(components):
                def_func = components
            if not components or callable(components):
                components_idx = list(range(0, array.shape[0]))
                for c_idx in components_idx:
                    components_dic[c_idx] = def_func
            for k,v in list(components_dic.items()):
                if v is None:
                    components_dic[k] = def_func
            #if self.overlap():
            if self.is_proc_inside(portion):
                ncx, ncy, ncz = self.get_no_cells(portion)
                if (ncx, ncy, ncz) != array.shape[1:]:
                    print ("self-Error in dump_region(): array and processor portion of different size.")
                    MPI.COMM_WORLD.Abort(errorcode=1)

                for i in range(portion[0], portion[1]+1):
                    for j in range(portion[2], portion[3]+1):
                        for k in range(portion[4], portion[5]+1):
                            cell_coords = self.map_cell2coord(i, j, k)
                            if coords != "mine":
                                if self.realm == CPL.CFD_REALM:
                                    cell_coords = self.map_cfd2md_coord(cell_coords)
                                else:
                                  cell_coords = self.map_md2cfd_coord(cell_coords)
                            [i_loc, j_loc, k_loc] = self.map_glob2loc_cell(portion, [i, j, k])
                            lines += str(cell_coords[0] + dx/2.0) + " "\
                                   + str(cell_coords[1] + dy/2.0) + " "\
                                   + str(cell_coords[2] + dz/2.0)
                                   
                            for k, f in list(components_dic.items()):
                                lines += " " + str(f(array[k, i_loc, j_loc, k_loc]))
                            lines += "\n"

            # Gather all the forces from every processor and dump them to a file at the root
            lines = comm.gather(lines, root=0)

            myrank = comm.Get_rank()
            if myrank == 0:
                with open(fname, "w") as file_out:
                    file_out.writelines(lines)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def cart_create(old_comm, dims, periods, coords):
    dummy_cart_comm = old_comm.Create_cart(dims, periods)
    temp_comm = old_comm.Split(0, dummy_cart_comm.Get_rank())
    new_cart_comm = temp_comm.Create_cart(dims, periods)
    comm_coords = new_cart_comm.Get_coords(new_cart_comm.Get_rank())
    if (not (coords == comm_coords).all()):
            print ("cart_create Error")
            exit()
    return new_cart_comm


# -----------------------------TESTING ROUTINES------------------------------ #


CONFIG_FILE = "COUPLER.in"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_NAME = os.path.basename(os.path.realpath(__file__))
TESTS_DIR_NAMES = ["initialisation", "mapping"]

def copyanything(src_dir, dst_dir, name):
    src_dir = os.path.join(src_dir, name)
    try:
        shutil.copytree(src_dir, os.path.join(dst_dir, name))
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src_dir, dst_dir)
        else: raise

def parametrize_file(source, dest, params):
    with open(source, "r+") as param_file_in:
        lines = param_file_in.readlines()
        for (k, v) in list(params.items()):
            lines = [l.replace("$[" + str(k) + "]", str(v)) for l in lines]
    with open(dest, "w+") as param_file_out:
        param_file_out.writelines(lines)

def parametrize_config(template_dir, params):
    # It assumes is in the temp directory with cpl/ folder accessible
    # from this level.
    source = os.path.join(template_dir, CONFIG_FILE)
    dest = os.path.join("cpl/", CONFIG_FILE)
    parametrize_file(source, dest, params)


def prepare_config(tmpdir, test_dir, md_fname, cfd_fname):
    tmpdir.mkdir("cpl")
    copyanything(test_dir, tmpdir.strpath, md_fname)
    copyanything(test_dir, tmpdir.strpath, cfd_fname)
    os.chdir(tmpdir.strpath)


def run_test(template_dir, config_params, md_exec, md_fname, md_args, cfd_exec,
             cfd_fname, cfd_args, md_params, cfd_params, err_msg, 
             debug=False, mpirun="split", printoutput=False):

    from distutils.spawn import find_executable
    import sys
    parametrize_config(template_dir, config_params)
    #Save parameter dictonaries to be read by md/cfd codes
    pickle.dump(md_params, open("md_params.dic", "wb"))
    pickle.dump(cfd_params, open("cfd_params.dic", "wb"))
    try:
        mdprocs = md_params["npx"] * md_params["npy"] * md_params["npz"]
        cfdprocs = cfd_params["npx"] * cfd_params["npy"] * cfd_params["npz"]
        if find_executable("mpiexec") is None:
            print("Error: mpiexec not found.")
            sys.exit(1)
        if find_executable(md_exec)is None:
            print(("Error: %s not found." % md_exec))
            sys.exit(1)
        if find_executable(cfd_exec)is None:
            print(("Error: %s not found." % cfd_exec))
            sys.exit(1)

            
        if os.path.exists(md_fname) and os.path.exists(cfd_fname):

            if "port" in mpirun:
                # For OpenMPI >= 3.0.0
                # cmd = " ".join(["mpiexec --oversubscribe", "-n", str(mdprocs), md_exec, md_args, 
                #                 "& PID=$!;", "mpiexec --oversubscribe", "-n", str(cfdprocs), 
                cmd = " ".join(["mpiexec", "-n", str(mdprocs), md_exec, md_args, 
                                "& PID=$!;", "mpiexec", "-n", str(cfdprocs), 
                                cfd_exec, cfd_args, "; wait $PID"])
#                cmd = " ".join(["cplexec",
#                                "-m ", str(mdprocs), " ' ", md_exec, md_args, " ' "
#                                "-c ", str(cfdprocs), " ' ", cfd_exec, cfd_args, " ' "])
            elif "split" in mpirun:
                # For OpenMPI >= 3.0.0
                # cmd = " ".join(["mpiexec --oversubscribe", "-n", str(mdprocs), md_exec, md_args,
                cmd = " ".join(["mpiexec","-n", str(mdprocs), md_exec, md_args,
                            ":", "-n", str(cfdprocs), cfd_exec, cfd_args])
            else:
                raise ValueError("MPIrun type unknown", mpirun)

            #Check for OpenMPI version greater than 3 and add oversubscribe option
            if (CPL.MPI_version == "OPENMPI"):
                if (CPL.ompi_major_version_no >= 3):
                    cmd.replace("mpiexec","mpiexec --oversubscribe")

            cmd.replace("mpiexec","mpiexec --oversubscribe")

            if debug:
                print(("\nMPI run: " + cmd))
            out = check_output(cmd, stderr=STDOUT, shell=True)
            if printoutput:
                print(out)
        else:
            print(("Current directory: " + os.getcwd()))
            print((md_fname + " or " + cfd_fname + " are not found."))
            assert False
            return False

    #This checsk the error message is as expected
    except CalledProcessError as exc:
        print((exc.output))
        if err_msg != "":
            print(("ERROR = ", err_msg))
            assert err_msg in exc.output.decode("utf-8")
        else:
            assert exc.output.decode("utf-8") == ""
    else:
        if err_msg != "":
            assert False
        else:
            assert True
    return True

def exec_tests(test="all"):
    import pytest
    import os
    test_path = os.path.dirname(os.path.realpath(__file__))
    test_path = os.path.join(test_path, "test")
    print(test_path)
    if test != "all":
        test_path = os.path.join(test_path, test)
    pytest.main(["-v", test_path])

def get_test_dir():
    import os
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test")

if __name__ == "__main__":
    lib = CPL()
