import pytest
import os
from subprocess import STDOUT, check_output, CalledProcessError
import shutil
import cPickle

TESTS_TEMPLATE_DIR = os.path.join(os.environ["CPLPY_PATH"], "tests/templates")
CONFIG_FILE = "COUPLER.in"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_NAME = os.path.basename(os.path.realpath(__file__))
MD_FILE = "md_test.py"
CFD_FILE = "cfd_test.py"


def parametrizeConfig(params):
    # It assumes is in the temp directory with cpl/ folder accessible
    # from this level.
    with open(os.path.join(TESTS_TEMPLATE_DIR,
                           CONFIG_FILE), "r+") as config_file:
        lines = config_file.readlines()
        for (k, v) in params.items():
            lines = [l.replace("$[" + str(k) + "]", str(v)) for l in lines]
    with open(os.path.join("cpl/", CONFIG_FILE), "w") as config_file:
        config_file.writelines(lines)



@pytest.fixture()
def prepareConfigEvery(tmpdir):
    tmpdir.mkdir("cpl")
    shutil.copy(os.path.join(TEST_DIR, MD_FILE), tmpdir.strpath)
    shutil.copy(os.path.join(TEST_DIR, CFD_FILE), tmpdir.strpath)
    os.chdir(tmpdir.strpath)
    os.chdir(tmpdir.strpath)
    return tmpdir


def run_test(config_params, md_params, cfd_params, err_msg):
    parametrizeConfig(config_params)
    cPickle.dump(md_params, open("md_params.dic", "wb"))
    cPickle.dump(cfd_params, open("cfd_params.dic", "wb"))
    try:
        mdprocs = md_params["npx"] * md_params["npy"] * md_params["npz"]
        cfdprocs = cfd_params["npx"] * cfd_params["npy"] * cfd_params["npz"]
        cmd = " ".join(["mpiexec", "-n", str(mdprocs), "python",
                        "md_test.py", ":", "-n", str(cfdprocs), "python",
                        "cfd_test.py"])

        check_output(cmd, stderr=STDOUT, shell=True)

    except CalledProcessError as exc:
        print exc.output
        if err_msg != "":
            assert err_msg in exc.output
        else:
            assert exc.output == ""
    else:
        if err_msg != "":
            assert False
        else:
            assert True

# -----INITIALISATION TESTS-----

# [TESTS WITH EXPECTED ERROR]

# EXPLANATION: These tests fail due to no_procs(MD) < no_procs(CFD) in
#              one direction.
@pytest.mark.parametrize("mdprocs,err_msg", [
                         ((0, 0, 0), "CFD or MD realm is missing"),
                         ((1, 4, 4), "number of MD processors in x must be"),
                         ((4, 4, 2), "number of MD processors in z must be"),
                         ((4, 3, 4), "number of MD processors in y must be"),
                         ((4, 4, 4), "")])
def test_mdprocs(prepareConfigEvery, mdprocs, err_msg):
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

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, err_msg)


# EXPLANATION: These tests fail due to no_procs(MD) != k*no_procs(CFD),
#              k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("cfdprocs, err_msg", [
                         ((2, 2, 3), "number of MD processors in z must be an integer"),
                         ((3, 2, 2), "number of MD processors in x must be an integer"),
                         ((2, 3, 2), ""),
                         ((4, 4, 6), "number of MD processors in z must be an integer"),
                         ((4, 6, 4), "number of MD processors in y must be greater"),
                         ((6, 4, 4), "number of MD processors in x must be an integer")])
def test_cfdprocs(prepareConfigEvery, cfdprocs, err_msg):
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

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, err_msg)


# EXPLANATION: These tests fail due to bad ranges in overlap cell ranges OR
#              cells in y-direction inside overlap region spanning across more
#              than one CFD processor.
@pytest.mark.parametrize("olapcells,expect_error", [
                         ((0, 0, 0, 0, 0, 0), "Overlap region has been specified outside"),
                         ((1, 25, 1, 12, 1, 24), "Overlap region is larger than the MD region"),
                         ((1, 24, 1, 12, 1, 25), "Overlap region is larger than the MD region"),
                         ((0, 24, 1, 12, 1, 24), "Overlap region is larger than the MD region"),
                         ((1, 24, 1, 12, 0, 24), "Overlap region is larger than the MD region"),
                         ((-2, 14, 1, 12, 1, 24), "Overlap region limits contains a negative index"),
                         ((1, 24, 1, 12, -2, 14), "Overlap region limits contains a negative index"),
                         ((1, 24, 1, -12, 1, 24), "Overlap region limits contains a negative index"),
                         ((15, 5, 1, 12, 1, 24), "Overlap region lower limits are greater than")])
def test_olapcells(prepareConfigEvery, olapcells, expect_error):
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

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)


# EXPLANATION: These tests fails due to k*no_cells != no_procs(CFD) OR,
#              k*no_cells != no_procs(MD), k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("domaincells,expect_error", [
                         ((1, 1, 1), "cells in the cfd domain is not an integer multipl"),
                         ((3, 4, 4), "cells in the cfd domain is not an integer multiple"),
                         ((4, 4, 3), "cells in the cfd domain is not an integer multiple"),
                         ((8, 8, 5), "cells in the cfd domain is not an integer multiple"),
                    #TODO     ((8, 3, 8), ""),
                         ((5, 8, 8), "cells in the cfd domain is not an integer multipl")])
def test_domaincells(prepareConfigEvery, domaincells, expect_error):
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

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)
