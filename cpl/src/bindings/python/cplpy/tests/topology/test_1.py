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
RUN_FILE = "run.sh"


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


def prepareConfig(tmp_path):
    shutil.copy(os.path.join(TEST_DIR, MD_FILE), tmp_path.strpath)
    shutil.copy(os.path.join(TEST_DIR, CFD_FILE), tmp_path.strpath)
    shutil.copy(os.path.join(TEST_DIR, RUN_FILE), tmp_path.strpath)
    os.chdir(tmp_path.strpath)


@pytest.fixture(scope="module")
def prepareConfigOnce(tmpdir_factory):
    tmpdir_factory.mktemp("cpl", numbered=False)
    prepareConfig(tmpdir_factory.getbasetemp())
    return tmpdir_factory.getbasetemp()


@pytest.fixture()
def prepareConfigEvery(tmpdir):
    tmpdir.mkdir("cpl")
    prepareConfig(tmpdir)
    os.chdir(tmpdir.strpath)
    return tmpdir


def run_test(config_params, md_params, cfd_params, err):
    parametrizeConfig(config_params)
    cPickle.dump(md_params, open("md_params.dic", "wb"))
    cPickle.dump(cfd_params, open("cfd_params.dic", "wb"))
    try:
        mdprocs = md_params["npx"] * md_params["npy"] * md_params["npz"]
        cfdprocs = cfd_params["npx"] * cfd_params["npy"] * cfd_params["npz"]
        check_output(["sh", RUN_FILE, str(mdprocs), str(cfdprocs)],
                     stderr=STDOUT)
    except CalledProcessError as exc:
        if err:
            assert exc.output != ""
        else:
            assert exc.output == ""
    else:
        if err:
            assert False
        else:
            assert True

# -----INITIALISATION TESTS-----

# [TESTS WITH EXPECTED ERROR]

# EXPLANATION: These tests fail due to no_procs(MD) < no_procs(CFD) in
#              one direction.
@pytest.mark.parametrize("mdprocs,expect_error", [
                         ((0, 0, 0), True),
                         ((1, 1, 1), True),
                         ((2, 2, 2), True),
                         ((3, 3, 3), True)])
def test_mdprocs(prepareConfigEvery, mdprocs, expect_error):
    MD_PARAMS = {"npx": 0, "npy": 0, "npz": 0,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "tstep_ratio": 50, }

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)


# EXPLANATION: These tests fail due to no_procs(MD) != k*no_procs(CFD),
#              k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("cfdprocs,expect_error", [
                         ((0, 0, 0), True),
                         ((3, 3, 3), True),
                         ((4, 4, 5), True),
                         ((4, 5, 4), True),
                         ((5, 4, 4), True)])
def test_cfdprocs(prepareConfigEvery, cfdprocs, expect_error):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"npx": 0, "npy": 0, "npz": 0,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "tstep_ratio": 50, }

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)


# EXPLANATION: These tests fail due to bad ranges in overlap cell ranges OR
#              cells in y-direction inside overlap region spanning across more
#              than one CFD processor.
@pytest.mark.parametrize("olapcells,expect_error", [
                         ((0, 0, 0, 0, 0, 0), True),
                         ((1, 25, 1, 12, 1, 24), True),
                         ((1, 24, 1, 12, 1, 25), True),
                         ((0, 24, 1, 12, 1, 24), True),
                         ((1, 24, 1, 12, 0, 24), True),
                         ((-2, 24, 1, 12, 1, 24), True),
                         ((1, 24, 1, 12, -2, 24), True),
                         ((1, 24, 1, -12, -2, 24), True)])
# This should fail       ((15, 5, 1, 12, 1, 24), True)])
def test_olapcells(prepareConfigEvery, olapcells, expect_error):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"npx": 2, "npy": 2, "npz": 2,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 0, "olap_xhi": 0,
                     "olap_ylo": 0, "olap_yhi": 0,
                     "olap_zlo": 0, "olap_zhi": 0,
                     "tstep_ratio": 50, }
    CONFIG_PARAMS["olap_xlo"], CONFIG_PARAMS["olap_xhi"] = olapcells[0:2]
    CONFIG_PARAMS["olap_ylo"], CONFIG_PARAMS["olap_yhi"] = olapcells[2:4]
    CONFIG_PARAMS["olap_zlo"], CONFIG_PARAMS["olap_zhi"] = olapcells[4:6]

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)


# EXPLANATION: These tests fails due to no_cells != k*no_procs(CFD) OR,
#              no_cells != k*no_procs(MD), k in [1,2,3,...] in one direction.
@pytest.mark.parametrize("domaincells,expect_error", [
                         ((0, 0, 0), True),
                         ((3, 3, 3), True),
                         ((4, 4, 5), True),
                         ((4, 5, 4), True),
                         ((5, 4, 4), True)])
def test_domaincells(prepareConfigEvery, domaincells, expect_error):
    MD_PARAMS = {"npx": 4, "npy": 4, "npz": 4,
                 "lx": 24.0, "ly": 24.0, "lz": 24.0, }

    CFD_PARAMS = {"npx": 0, "npy": 0, "npz": 0,
                  "lx": 24.0, "ly": 24.0, "lz": 24.0,
                  "ncx": 24, "ncy": 24, "ncz": 24, }
    CFD_PARAMS["ncx"], CFD_PARAMS["ncy"], CFD_PARAMS["ncz"] = domaincells

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 24,
                     "olap_ylo": 1, "olap_yhi": 4,
                     "olap_zlo": 1, "olap_zhi": 24,
                     "tstep_ratio": 50, }

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)
