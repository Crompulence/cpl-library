#!/usr/bin/env python2
import pytest
import cplpy
from cplpy import run_test, prepare_config, get_test_dir
import os
import itertools
import numpy as np

# -----MAPPING TESTS-----

# EXPLANATION: These tests fail due to no_procs(MD) != k*no_procs(CFD),
#              k in [1,2,3,...] in one direction.

MD_FNAME = "md_test.py"
MD_ARGS = MD_FNAME
MD_EXEC = "python2"
CFD_FNAME = "cfd_test.py"
CFD_ARGS = CFD_FNAME
CFD_EXEC = "python2"
TEST_TEMPLATE_DIR = os.path.join(get_test_dir(), "templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)

#Get all permutations
ncx = 24; ncy = 24; ncz = 24
maxprocperdir = 6
maxprocs = 16
tests = ["cell_test", "coord_test"]
if (cplpy.CPL.MPI_version == "OPENMPI"):
    cases=["split"]
else:
    cases=["port", "split"]
cases = list(itertools.combinations_with_replacement(range(1,maxprocperdir), 3))
perms = []; n=0
for cpltype in cpltypes:
    for test in tests:
        for i in range(len(cases)):
            for j in range(len(cases)-i):
                if ((cases[i+j][0] < cases[i][0]) or
                    (cases[i+j][1] < cases[i][1]) or
                    (cases[i+j][2] < cases[i][2])):
                    errstr="number of MD processors must be greater than or equal to CFD processors"
                elif ((ncx%cases[i+j][0] != 0) or 
                      (ncy%cases[i+j][1] != 0) or
                      (ncz%cases[i+j][2] != 0)):
                    errstr="The number of cells in the cfd domain is not an integer multiple of the number of processors"
                elif ((cases[i+j][0]%cases[i][0] != 0) or
                      (cases[i+j][2]%cases[i][2] != 0)):
                    errstr="number of MD processors must be an integer multiple of number of CFD processors"
                else:
                    errstr=""
                #Skip any errors greater than 16 processors and cases greater than 32
                if (np.product(cases[i])+np.product(cases[i+j]) < maxprocs):
                    perms.append((cpltype, test, cases[i], cases[i+j], errstr))
                    print(n, perms[-1])
                    n += 1



@pytest.mark.parametrize("cpltype, test, cfdprocs, mdprocs, err_msg", perms)
def test_mapcells(prepare_config_fix, cpltype, test, cfdprocs, mdprocs, err_msg):
    MD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
                 "which_test": test}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": ncx, "ncy": ncy, "ncz": ncz,
                  "which_test": "cell_test"}
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "cnst_xlo": 1, "cnst_xhi": 1,
                     "cnst_ylo": 1, "cnst_yhi": 1,
                     "cnst_zlo": 1, "cnst_zhi": 1,
                     "tstep_ratio": 50, }

    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
             CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, 
             mpirun=cpltype)


#@pytest.mark.parametrize("cfdprocs, mdprocs, err_msg", perms)
#def test_mappoint(prepare_config_fix, cfdprocs, mdprocs, err_msg):
#    MD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
#                 "which_test": test}
#    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

#    CFD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
#                  "ncx": ncx, "ncy": ncy, "ncz": ncz,
#                  "which_test": "coord_test"}
#    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

#    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
#                     "olap_xlo": 1, "olap_xhi": 24,
#                     "olap_ylo": 1, "olap_yhi": 4,
#                     "olap_zlo": 1, "olap_zhi": 24,
#                     "cnst_xlo": 1, "cnst_xhi": 1,
#                     "cnst_ylo": 1, "cnst_yhi": 1,
#                     "cnst_zlo": 1, "cnst_zhi": 1,
#                     "tstep_ratio": 50, }

#    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
#             CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg)
