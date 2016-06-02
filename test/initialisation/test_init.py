import pytest
from cplpy import run_test, prepare_config
import os


# -----INITIALISATION TESTS-----

# [TESTS WITH EXPECTED ERROR]

# EXPLANATION: These tests fail due to no_procs(MD) < no_procs(CFD) in
#              one direction.

MD_FNAME = "md_test.py"
MD_RUN = "python " + MD_FNAME
CFD_FNAME = "cfd_test.py"
CFD_RUN = "python " + CFD_FNAME
TEST_TEMPLATE_DIR = os.path.join(os.environ["CPL_PATH"], "test/templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)


@pytest.mark.parametrize("mdprocs, err_msg", [
                         ((0, 0, 0), "CFD or MD realm is missing"),
                         ((1, 4, 4), "number of MD processors in x must be"),
                         ((4, 4, 2), "number of MD processors in z must be"),
                         ((4, 3, 4), "number of MD processors in y must be"),
                         ((4, 4, 4), "")])
def test_mdprocs(prepare_config_fix, mdprocs, err_msg):
    MD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "cnst_xlo": 1, "cnst_xhi": 1,
                     "cnst_ylo": 1, "cnst_yhi": 1,
                     "cnst_zlo": 1, "cnst_zhi": 1,
                     "tstep_ratio": 50, }

    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_RUN, CFD_RUN, MD_PARAMS,
             CFD_PARAMS, err_msg)


# EXPLANATION: These tests fail due to no_procs(MD) != k*no_procs(CFD),
#              k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("cfdprocs, err_msg", [
                         ((2, 2, 3), "number of MD processors in z must be an integer"),
                         ((3, 2, 2), "number of MD processors in x must be an integer"),
                         ((2, 3, 2), ""),
                         ((4, 4, 6), "number of MD processors in z must be an integer"),
                         ((4, 6, 4), "number of MD processors in y must be greater"),
                         ((6, 4, 4), "number of MD processors in x must be an integer")])
def test_cfdprocs(prepare_config_fix, cfdprocs, err_msg):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "cnst_xlo": 1, "cnst_xhi": 1,
                     "cnst_ylo": 1, "cnst_yhi": 1,
                     "cnst_zlo": 1, "cnst_zhi": 1,
                     "tstep_ratio": 50, }

    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_RUN, CFD_RUN, MD_PARAMS,
             CFD_PARAMS, err_msg)


# EXPLANATION: These tests fail due to bad ranges in overlap cell ranges OR
#              cells in y-direction inside overlap region spanning across more
#              than one CFD processor.
@pytest.mark.parametrize("olapcells, err_msg", [
                         ((0, 0, 0, 0, 0, 0), "Overlap region has been specified outside"),
                         ((1, 25, 1, 12, 1, 24), "Overlap region is larger than the MD region"),
                         ((1, 24, 1, 12, 1, 25), "Overlap region is larger than the MD region"),
                         ((0, 24, 1, 12, 1, 24), "Overlap region is larger than the MD region"),
                         ((1, 24, 1, 12, 0, 24), "Overlap region is larger than the MD region"),
                         ((-2, 14, 1, 12, 1, 24), "Overlap region limits contains a negative index"),
                         ((1, 24, 1, 12, -2, 14), "Overlap region limits contains a negative index"),
                         #TODO:Test fails with: Overlap region lower limits are greater than upper limits for some directions
                         #((1, 24, 1, -12, 1, 24), "Overlap region limits contains a negative index"),
                         ((15, 5, 1, 12, 1, 24), "Overlap region lower limits are greater than")])
def test_olapcells(prepare_config_fix, olapcells, err_msg):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"npx": 2, "npy": 2, "npz": 2,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "cnst_xlo": 1, "cnst_xhi": 1,
                     "cnst_ylo": 1, "cnst_yhi": 1,
                     "cnst_zlo": 1, "cnst_zhi": 1,
                     "tstep_ratio": 50, }
    CONFIG_PARAMS["olap_xlo"], CONFIG_PARAMS["olap_xhi"] = olapcells[0:2]
    CONFIG_PARAMS["olap_ylo"], CONFIG_PARAMS["olap_yhi"] = olapcells[2:4]
    CONFIG_PARAMS["olap_zlo"], CONFIG_PARAMS["olap_zhi"] = olapcells[4:6]

    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_RUN, CFD_RUN, MD_PARAMS,
             CFD_PARAMS, err_msg)


# EXPLANATION: These tests fails due to k*no_cells != no_procs(CFD) OR,
#              k*no_cells != no_procs(MD), k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("domaincells, err_msg", [
                         ((1, 1, 1), "cells in the cfd domain is not an integer multipl"),
                         ((3, 4, 4), "cells in the cfd domain is not an integer multiple"),
                         ((4, 4, 3), "cells in the cfd domain is not an integer multiple"),
                         ((8, 8, 5), "cells in the cfd domain is not an integer multiple"),
                    #TODO     ((8, 3, 8), ""),
                         ((5, 8, 8), "cells in the cfd domain is not an integer multipl")])
def test_domaincells(prepare_config_fix, domaincells, err_msg):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"npx": 1, "npy": 1, "npz": 1,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0}
    CFD_PARAMS["ncx"], CFD_PARAMS["ncy"], CFD_PARAMS["ncz"] = domaincells

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 1,
                     "olap_ylo": 1, "olap_yhi": 1,
                     "olap_zlo": 1, "olap_zhi": 1,
                     "cnst_xlo": 1, "cnst_xhi": 1,
                     "cnst_ylo": 1, "cnst_yhi": 1,
                     "cnst_zlo": 1, "cnst_zhi": 1,
                     "tstep_ratio": 50, }

    run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_RUN, CFD_RUN, MD_PARAMS,
             CFD_PARAMS, err_msg)
