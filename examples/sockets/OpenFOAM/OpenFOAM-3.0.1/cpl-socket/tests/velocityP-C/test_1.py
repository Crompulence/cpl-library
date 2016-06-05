import pytest
import os
from subprocess import STDOUT, check_output, CalledProcessError
import shutil
import cPickle
import numpy as np

# Template for the coupler in put file from cplpy tests/templates folder
TESTS_TEMPLATE_DIR = os.path.join(os.environ["CPLPY_PATH"], "tests/templates")
CONFIG_FILE = "COUPLER.in"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_NAME = os.path.basename(os.path.realpath(__file__))
MD_FILE = "lammps/vels.in"
CFD_FILE = "cfd_vels.py"


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


def prepare(tmp_path):
    shutil.copy(os.path.join(TEST_DIR, MD_FILE), tmp_path.strpath)
    shutil.copy(os.path.join(TEST_DIR, CFD_FILE), tmp_path.strpath)
    os.chdir(tmp_path.strpath)


@pytest.fixture()
def prepareConfig(tmpdir):
    tmpdir.mkdir("cpl")
    prepare(tmpdir)
    os.chdir(tmpdir.strpath)
    return tmpdir


#TODO: Add error string to check for it in error output
def run_test(config_params, md_params, cfd_params, err):
    parametrizeConfig(config_params)
    cPickle.dump(cfd_params, open("cfd_params.dic", "wb"))
    try:
        mdprocs = md_params["npx"] * md_params["npy"] * md_params["npz"]
        cfdprocs = cfd_params["npx"] * cfd_params["npy"] * cfd_params["npz"]
        check_output(" ".join(["mpiexec", "-n", str(mdprocs), "lmp_cpl",
                     "vels.in", ":", "-n", str(cfdprocs), "python",
                     "cfd_vels.py"]), stderr=STDOUT, shell=True)
    except CalledProcessError as exc:
        if err:
            assert exc.output != ""
        else:
            assert exc.output == ""
    else:
        if err:
            assert False



def compare_vels(tol, lammps_fname="lammps_vels.dat", cfd_fname="cfd_vels.dat"):

    # Line format of CFD script file -- > x z vx vy vz
    with open(cfd_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
    cfd_cells = {}
    for l in cfd_lines:
        cfd_cells[(float(l[0]), float(l[1]))] = np.array([float(l[2]),
                                                         float(l[3]),
                                                         float(l[4])])

    # Line format of LAMMPS file -- > chunk x z ncount vx vy vz
    with open(lammps_fname) as lammps_file:
        lammps_lines = lammps_file.readlines()
    skip = int(lammps_lines[3].split(" ")[1])
    lammps_lines = lammps_lines[4:]
    lammps_lines = lammps_lines[:skip]
    lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
    lammps_cells = {}
    for l in lammps_lines:
        l = filter(None, l)
        lammps_cells[(float(l[1]), float(l[2]))] = np.array([float(l[4]),
                                                            float(l[5]),
                                                            float(l[6])])
    # Check if the number of cells are the same 
    if len(lammps_cells) != len(cfd_cells):
        assert False

    # Compare each cell velocity up to a certain tolerance
    for cell in lammps_cells.keys():
        try:
            diff_vel = abs(cfd_cells[cell] - lammps_cells[cell])
            if (np.any(diff_vel > tol)):
                print "False"
                print cfd_cells[cell]
                print lammps_cells[cell]
                assert False
        except KeyError:
            print "Cell not found: "
            print cell
            assert False

# -----VELOCITY TESTS-----

# EXPLANATION: These tests fail due to no_procs(MD) < no_procs(CFD) in
#              one direction.
@pytest.mark.parametrize("cfdprocs, mdprocs, expect_error", [
                         ((3, 3, 3), (3, 3, 3),  False)])
#                         ((5, 4, 4), (5, 4, 4),  True)])
def test_velocitiesP2C(prepareConfig, cfdprocs, mdprocs, expect_error):
    MD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0,
                  "ncx": 15, "ncy": 15, "ncz": 15, }
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 15,
                     "olap_ylo": 1, "olap_yhi": 5,
                     "olap_zlo": 1, "olap_zhi": 15,
                     "cnst_xlo": 1, "cnst_xhi": 15,
                     "cnst_ylo": 5, "cnst_yhi": 5,
                     "cnst_zlo": 1, "cnst_zhi": 15,
                     "tstep_ratio": 1, }

    run_test(CONFIG_PARAMS, MD_PARAMS, CFD_PARAMS, expect_error)
    compare_vels(1e-6)


