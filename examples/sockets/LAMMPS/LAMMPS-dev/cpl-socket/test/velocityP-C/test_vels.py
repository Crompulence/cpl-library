import pytest
import numpy as np
from cplpy import run_test, prepare_config
import os

# -----Velocities TESTS-----

# EXPLANATION:

MD_FNAME = "lammps_vels.in"
MD_ARGS = "-in " + MD_FNAME
MD_EXEC = "lmp_cpl"
CFD_FNAME = "dummyCFD.py"
CFD_ARGS = CFD_FNAME
CFD_EXEC = "python"
TEST_TEMPLATE_DIR = os.path.join(os.environ["CPL_PATH"], "test/templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)


def compare_vels(tol, lammps_fname="lammps_vels.dat",
                 cfd_fname="cfd_vels.dat"):

    # Line format of CFD script file -- > x y z vx vy vz
    with open(cfd_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
    cfd_cells = {}
    for l in cfd_lines:
        cfd_cells[(float(l[0]), float(l[1]), float(l[2]))] = np.array([float(l[3]),
                                                        	       float(l[4]),
                                                         	       float(l[5])])

    # Line format of LAMMPS file -- > chunk x y z ncount vx vy vz
    with open(lammps_fname, "r") as lammps_file:
        lammps_lines = lammps_file.readlines()
    skip = int(lammps_lines[3].split(" ")[1])
    lammps_lines = lammps_lines[4:]
    lammps_lines = lammps_lines[:skip]
    lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
    lammps_cells = {}
    for l in lammps_lines:
        l = filter(None, l)
        lammps_cells[(float(l[1]), float(l[2]), float(l[3]))] = np.array([float(l[5]),
                                                            		  float(l[6]),
                                                            		  float(l[7])])

    # Compare each cell velocity up to a certain tolerance
    for cell in cfd_cells.keys():
        try:
            diff_vel = abs(cfd_cells[cell] - lammps_cells[cell])
            if (np.any(diff_vel > tol)):
                print "Cell value disagreement:"
                print cfd_cells[cell]
                print lammps_cells[cell]
                assert False
        except KeyError:
            print "Cell not found: " + str(cell)
            assert False

# -----VELOCITY TESTS-----

# EXPLANATION: See README-test located in this folder.


@pytest.mark.parametrize("cfdprocs, mdprocs, err_msg", [
                         ((3, 3, 3), (3, 3, 3),  ""),
                         ((1, 1, 1), (3, 3, 3),  "")])
def test_velocitiesP2C(prepare_config_fix, cfdprocs, mdprocs, err_msg):
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

    correct = run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
                       CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, True)
    if correct:
        compare_vels(1e-6)
