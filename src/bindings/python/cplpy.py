from __future__ import print_function, division
from ctypes import c_char_p, c_char, c_int, c_double, c_bool, c_void_p, byref, POINTER, util, pointer
import ctypes
from mpi4py import MPI
import numpy as np
from numpy.ctypeslib import ndpointer, load_library
from functools import wraps
import os
import time
from subprocess import STDOUT, check_output, CalledProcessError
import shutil
import cPickle
import errno


__all__ = ["CPL", "cart_create", "run_test", "prepare_config", "parametrize_file"]


class OpenMPI_Not_Supported(Exception):
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
                 "ncx": c_int, "ncy": c_int, "ncz": c_int,
                 "npx_md": c_int, "npy_md": c_int, "npz_md": c_int,
                 "npx_cfd": c_int, "npy_cfd": c_int, "npz_cfd": c_int,
                 "overlap": c_int, "xl_md": c_double, "yl_md": c_double,
                 "nsteps_md": c_int, "nsteps_cfd": c_int, "nsteps_coupled": c_int,
                 "zl_md": c_double, "xl_cfd": c_double, "yl_cfd": c_double,
                 "zl_cfd": c_double, "dx" : c_double, "dy" : c_double, "dz" : c_double,
                 "x_orig_cfd": c_double,"y_orig_cfd": c_double,"z_orig_cfd": c_double,
                 "x_orig_md": c_double,"y_orig_md": c_double,"z_orig_md": c_double
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
    _lib_path = os.environ.get("CPL_LIBRARY_PATH")
    # Try using CPL_LIBRARY_PATH and if not look
    # in system path (or ctype path)
    try:
        _cpl_lib = load_library(_libname, _lib_path)
    except OSError as e:
        print("OSError: ", e)
        print("WARNING - " + _libname + " not found under path specified by \
              CPL_LIBRARY_PATH variable: " + _lib_path)
        print("Attempting to find system version")
        _lib_path = util.find_library(_libname)
    # If not found, try to look around the current directory system
    # for the library
    else:
        try:
            _cpl_lib = load_library(_libname, _lib_path)
        except OSError:
            print("WARNING - " + _libname + " not found in system path or \
                  $CPL_LIBRARY_PATH: " + _lib_path)
            print("Attempting to search current folder system")
            trydir = ''
            found = False
            for level in range(10):
                if not found:
                    trydir += "../"
                for root, dirs, files in os.walk(trydir):
                    print(level, trydir, root)
                    if "libcpl.so" in files:
                        print("A version of " + _libname + " is found under \
                              path: " + root + ". Attempting to use this..")
                        found = True
                        _lib_path = root
                        break

            if found:
                _cpl_lib = load_library(_libname, _lib_path)
            else:
                print(_libname + " not found, please ensure you have \
                      compiled correctly")
                print(" and set CPL_LIBRARY_PATH appropitely")
                time.sleep(2)
                MPI.COMM_WORLD.Abort(errorcode=1)

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

    _py_init = _cpl_lib.CPLC_init

    #OpenMPI comm greater than c_int
    if MPI._sizeof(MPI.Comm) == ctypes.sizeof(c_int): 
        _py_init.argtypes = [c_int, POINTER(c_int)]
    else:
        excptstr ="Problem is in create_comm wrapper, as the OpenMPI COMM handle is not "
        excptstr += "an integer, c_void_p should be used so C bindings needs something like **void" 
        excptstr += "(No idea what to do in the Fortran code, maybe MPI_COMM_f2C required)"
        raise OpenMPI_Not_Supported(excptstr)
        _py_init.argtypes = [c_int, POINTER(c_void_p)]

    @abortMPI
    def init(self, calling_realm):

        self.realm = calling_realm
        # Build a communicator mpi4py python object from the
        # handle returned by the CPL_init function.
        if MPI._sizeof(MPI.Comm) == ctypes.sizeof(c_int):
            MPI_Comm = c_int
        else:
            MPI_Comm = c_void_p

        # Call create comm
        returned_realm_comm = c_int()
        self._py_init(calling_realm, byref(returned_realm_comm))

        # Use an intracomm object as the template and override value
        newcomm = MPI.Intracomm()
        newcomm_ptr = MPI._addressof(newcomm)
        comm_val = MPI_Comm.from_address(newcomm_ptr)
        comm_val.value = returned_realm_comm.value

        return newcomm


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
            print ("CPL-ERROR: CPL Library function '" +
                   str(fun_name) + "' not found!")
            raise KeyError
        else:
            self._var = var_ctype()

            if ("array" in fun_name):
                print ("ENTRO")
                var_len = c_int()
                fun.argtypes.append(POINTER(c_int))
                print ("EY")
                fun(c_char_p(section), c_char_p(var_name), byref(self._var), byref(var_len))
                print ("len:" , var_len.value)
                #print (self._var[0])
                #print (byref(var[0]))
                a = ([self._var[i] for i in xrange(var_len.value)])
                return a
            else:
                fun(c_char_p(section), c_char_p(var_name), byref(self._var))
                return self._var.value



    _py_finalize = _cpl_lib.CPLC_finalize

    @abortMPI
    def finalize(self):
        self._py_finalize()

    py_setup_cfd = _cpl_lib.CPLC_setup_cfd

    py_setup_cfd.argtypes = \
        [c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def setup_cfd(self, icomm_grid, xyzL, 
                        xyz_orig, ncxyz):
        """
        Keyword arguments:
        real -- the real part (default 0.0)
        imag -- the imaginary part (default 0.0)
        """

        self.py_setup_cfd(MPI._handleof(icomm_grid), xyzL,
                          xyz_orig, ncxyz)


    py_setup_md = _cpl_lib.CPLC_setup_md

    py_setup_md.argtypes = \
        [c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous'),
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def setup_md(self, icomm_grid, xyzL, xyz_orig):
        """
        setup_md(self, dt, icomm_grid, xyzL, xyz_orig)
        Keyword arguments:
        real -- the real part (default 0.0)
        imag -- the imaginary part (default 0.0)
        """
        self.py_setup_md(MPI._handleof(icomm_grid), xyzL, xyz_orig)

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

    py_my_proc_extents = _cpl_lib.CPLC_my_proc_extents
    py_my_proc_extents.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def my_proc_extents(self):
        extents = np.zeros(6, order='F', dtype=np.int32)
        self.py_my_proc_extents(extents)
        return extents

    py_proc_portion = _cpl_lib.CPLC_proc_portion
    py_proc_portion.argtypes = \
        [ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous'),
         c_int,
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def proc_portion(self, coord, realm, limits):
        coord = self._type_check(coord)
        limits = self._type_check(limits)
        portion = np.zeros(6, order='F', dtype=np.int32)
        self.py_proc_portion(coord, realm, limits, portion)
        return portion

    py_my_proc_portion = _cpl_lib.CPLC_my_proc_portion
    py_my_proc_portion.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous'),
         ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def my_proc_portion(self, limits):
        limits = self._type_check(limits)
        portion = np.zeros(6, order='F', dtype=np.int32)
        self.py_my_proc_portion(limits, portion)
        return portion

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

    py_map_cell2coord = _cpl_lib.CPLC_map_cell2coord
    py_map_cell2coord.argtypes = \
        [c_int, c_int, c_int,
         ndpointer(np.float64, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_cell2coord(self, i, j, k):
        coord = np.zeros(3, order='F', dtype=np.float64)
        self.py_map_cell2coord(i, j, k, coord)
        return coord

    py_map_coord2cell = _cpl_lib.CPLC_map_coord2cell
    py_map_coord2cell.argtypes = \
        [c_double, c_double, c_double,
         ndpointer(np.int32, shape=(3,), flags='aligned, f_contiguous')]

    @abortMPI
    def map_coord2cell(self, x, y, z):
        cell = np.zeros(3, order='F', dtype=np.int32)
        self.py_map_coord2cell(x, y, z, cell)
        return cell

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

    py_get_olap_limits = _cpl_lib.CPLC_get_olap_limits
    py_get_olap_limits.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_olap_limits(self):
        limits = np.zeros(6, order='F', dtype=np.int32)
        self.py_get_olap_limits(limits)
        return limits

    py_get_cnst_limits = _cpl_lib.CPLC_get_cnst_limits
    py_get_cnst_limits.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def get_cnst_limits(self):
        limits = np.zeros(6, order='F', dtype=np.int32)
        self.py_get_cnst_limits(limits)
        return limits

    py_set_timing = _cpl_lib.CPLC_set_timing
    py_set_timing.argtypes = \
        [c_int, c_int, c_double]

    @abortMPI
    def set_timing(self, initialstep, nsteps, dt):
        self.py_set_timing(initialstep, nsteps, dt)

    py_send = _cpl_lib.CPLC_send
    py_send.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'), 
         POINTER(c_bool)]

    @abortMPI
    def send(self, asend, limits):
        asend = self._type_check(asend)
        asend_shape = np.array(asend.shape, order='F', dtype=np.int32)
        send_flag = c_bool()
        self.py_send(asend, asend_shape, limits, byref(send_flag))
        return send_flag.value

    py_recv = _cpl_lib.CPLC_recv
    py_recv.argtypes = \
        [ndpointer(np.float64, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'),
         ndpointer(np.int32, ndim=1, flags='aligned, f_contiguous'), 
         POINTER(c_bool)]

    @abortMPI
    def recv(self, arecv, limits):
        arecv = self._type_check(arecv)
        arecv_shape = np.array(arecv.shape, order='F', dtype=np.int32)
        recv_flag = c_bool()
        self.py_recv(arecv, arecv_shape, limits, byref(recv_flag))
        return arecv, recv_flag.value

    py_overlap = _cpl_lib.CPLC_overlap
    py_overlap.argtypes = []

    @abortMPI
    def overlap(self):
        self.py_overlap.restype = c_bool
        return self.py_overlap()

    py_is_proc_inside = _cpl_lib.CPLC_is_proc_inside
    py_is_proc_inside.argtypes = \
        [ndpointer(np.int32, shape=(6,), flags='aligned, f_contiguous')]

    @abortMPI
    def is_proc_inside(self, region):
        self.py_is_proc_inside.restype = c_bool
        region = self._type_check(region)
        return self.py_is_proc_inside(region)

    @abortMPI
    def get(self, var_name):
        try:
            var_type = _CPL_GET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_" + var_name)
        except KeyError:
            print ("CPL-ERROR: CPL Library function '" +
                   str(var_name) + "' not found!")
            raise KeyError
        else:
            fun.restype = var_type
            return fun()

    @abortMPI
    def set(self, var_name, value):
        try:
            var_type = _CPL_SET_VARS[var_name]
            fun = getattr(self._cpl_lib, "CPLC_set_" + var_name)
        except KeyError:
            print ("CPL-ERROR: CPL Library function '" +
                   str(var_name) + "' not found!")
            raise KeyError
        else:
            fun.argtypes = [var_type]
            return fun(var_type(value))

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
                components_idx = list(xrange(0, array.shape[0]))
                for c_idx in components_idx:
                    components_dic[c_idx] = def_func
            for k,v in components_dic.items():
                if v is None:
                    components_dic[k] = def_func
            #if self.overlap():
            if self.is_proc_inside(portion):
                ncx, ncy, ncz = self.get_no_cells(portion)
                if (ncx, ncy, ncz) != array.shape[1:]:
                    print ("self-Error in dump_region(): array and processor portion of different size.")
                    MPI.COMM_WORLD.Abort(errorcode=1)

                for i in xrange(portion[0], portion[1]+1):
                    for j in xrange(portion[2], portion[3]+1):
                        for k in xrange(portion[4], portion[5]+1):
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
                                   
                            for k, f in components_dic.items():
                                lines += " " + str(f(array[k, i_loc, j_loc, k_loc]))
                            lines += "\n"

            # Gather all the forces from every processor and dump them to a file at the root
            lines = comm.gather(lines, root=0)

            myrank = comm.Get_rank()
            if myrank == 0:
                with open(fname, "w") as file_out:
                    file_out.writelines(lines)




def cart_create(old_comm, dims, periods, coords):
    dummy_cart_comm = old_comm.Create_cart(dims, periods)
    temp_comm = old_comm.Split(0, dummy_cart_comm.Get_rank())
    new_cart_comm = temp_comm.Create_cart(dims, periods)
    comm_coords = new_cart_comm.Get_coords(new_cart_comm.Get_rank())
    if (not (coords == comm_coords).all()):
            print ("Not good!")
            exit()
    return new_cart_comm


# -----------------------------TESTING ROUTINES------------------------------ #


CONFIG_FILE = "COUPLER.in"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_NAME = os.path.basename(os.path.realpath(__file__))

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
        for (k, v) in params.items():
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
             cfd_fname, cfd_args, md_params, cfd_params, err_msg, debug=False, mpirun="port"):
    parametrize_config(template_dir, config_params)
    cPickle.dump(md_params, open("md_params.dic", "wb"))
    cPickle.dump(cfd_params, open("cfd_params.dic", "wb"))
    try:
        mdprocs = md_params["npx"] * md_params["npy"] * md_params["npz"]
        cfdprocs = cfd_params["npx"] * cfd_params["npy"] * cfd_params["npz"]
        if os.path.exists(md_fname) and os.path.exists(cfd_fname):
            if mpirun == "port":
                cmd = " ".join(["mpiexec", "-n", str(mdprocs), md_exec, md_args, 
                                "& PID=$!;", "mpiexec", "-n", str(cfdprocs), 
                                cfd_exec, cfd_args, "; wait $PID"])
            else:
                cmd = " ".join(["mpiexec", "-n", str(mdprocs), md_exec, md_args,
                                ":", "-n", str(cfdprocs), cfd_exec, cfd_fname])
            if debug:
                print ("\nMPI run: " + cmd)
            check_output(cmd, stderr=STDOUT, shell=True)
        else:
            print ("Current directory: " + os.getcwd())
            print (md_fname + " or " + cfd_fname + " are not found.")
            assert False
            return False

    except CalledProcessError as exc:
        print (exc.output)
        if err_msg != "":
            assert err_msg in exc.output
        else:
            assert exc.output == ""
    else:
        if err_msg != "":
            assert False
        else:
            assert True
    return True

if __name__ == "__main__":
    lib = CPL()
