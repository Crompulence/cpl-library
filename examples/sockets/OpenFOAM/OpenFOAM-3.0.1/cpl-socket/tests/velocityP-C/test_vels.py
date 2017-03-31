import pytest
import numpy as np
from cplpy import run_test, prepare_config, parametrize_file
import os
import sys
import subprocess

try:
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
except:
    "Error: PyFoam package is required to run the tests"
    sys.exit()

# -----Velocity TESTS-----

# EXPLANATION:

MD_FNAME = "dummyMD_vels.py"
MD_ARGS = MD_FNAME
MD_EXEC = "python"
CFD_FNAME = "test_vels_case"
CFD_ARGS = "-parallel -case " + CFD_FNAME
CFD_EXEC = "CPLIcoFoam"
TEST_TEMPLATE_DIR = os.path.join(os.environ["CPL_PATH"], "test/templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)



def compare_forces(tol, cfd_params, md_fname="md_vels.dat",
                   openfoam_dir="test_vels_case"):

    vel = ParsedParameterFile(openfoam_dir + "/2/U")["boundaryField"]["CPLReceiveMD"]["value"]
    vx,vy,vz = zip(*[tuple(v) for v in vel])
    boundary = ParsedParameterFile(openfoam_dir + "/constant/polyMesh/boundary", 
                                   boundaryDict=True)
    boundary = dict(boundary[i:i+2] for i in range(0, len(boundary), 2))
    nFaces = boundary["CPLReceiveMD"]["nFaces"]
    startFace = boundary["CPLReceiveMD"]["startFace"]
    owner = ParsedParameterFile(openfoam_dir + "/constant/polyMesh/owner",
                                listDictWithHeader=True)
    cell_cx = ParsedParameterFile(openfoam_dir + "/2/ccx")["internalField"]
    cell_cy = ParsedParameterFile(openfoam_dir + "/2/ccy")["internalField"]
    cell_cz = ParsedParameterFile(openfoam_dir + "/2/ccz")["internalField"]

    openfoam_cells = {}
    cell_no = 0
    for cell_no in xrange(0, nFaces): 
        cell_owner = owner[cell_no+startFace]
        cell_coord = (float(cell_cx[cell_owner]), float(cell_cy[cell_owner]), 
                      float(cell_cz[cell_owner]))
        # Openfoam output cell centres with 6 decimal figures
        k = "{0:.5f}".format(cell_coord[0]), "{0:.5f}".format(cell_coord[1]),\
            "{0:.5f}".format(cell_coord[2])
        openfoam_cells[k] = np.array([float(vx[cell_no]), 
                                      float(vy[cell_no]), 
                                      float(vz[cell_no])])


    # Line format of dummy md forces file -- > x y z sxy syy szy
    with open(md_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    md_lines = [l[:-1].split(" ") for l in cfd_lines]
    md_cells = {}
    for l in md_lines:
        k = "{0:.5f}".format(float(l[0])), "{0:.5f}".format(float(l[1])), "{0:.5f}".format(float(l[2]))
        md_cells[k] = np.array([float(l[3]), float(l[4]), float(l[5])])

    
    for k in md_cells.keys():
        try:
            diff_forces = abs(md_cells[k] - openfoam_cells[k])
            if (np.any(diff_forces > tol)):
                print md_cells[k]
                print openfoam_cells[k]
                assert False
        except KeyError:
            print "Cell not found: cell " + k
            assert False



# -----FORCES TESTS-----

# EXPLANATION: See README-test located in this folder.

@pytest.mark.parametrize("cfdprocs, mdprocs, cells, err_msg", [
                         ((3, 3, 3), (3, 3, 3), (15, 15, 15), ""),
                         ((1, 2, 1), (3, 2, 2), (30, 36, 24), ""),
                         ((4, 3, 3), (4, 3, 3), (20, 15, 27), ""),
                         ((3, 3, 3), (3, 3, 3), (30, 15, 21), "")])
def test_velsP2C(prepare_config_fix, cfdprocs, mdprocs, cells,  err_msg):
    MD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs
    CFD_PARAMS["ncx"], CFD_PARAMS["ncy"], CFD_PARAMS["ncz"] = cells
    # Needed for decomposParDict
    CFD_PARAMS["nprocs"] = cfdprocs[0]*cfdprocs[1]*cfdprocs[2]

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": cells[0],
                     "olap_ylo": 1, "olap_yhi": 5,
                     "olap_zlo": 1, "olap_zhi": cells[2],
                     "cnst_xlo": 1, "cnst_xhi": cells[0],
                     "cnst_ylo": 3, "cnst_yhi": 5,
                     "cnst_zlo": 1, "cnst_zhi": cells[2],
                     "tstep_ratio": 1, }

    # Parametrize OpenFOAM files
    mesh_file = os.path.join(CFD_FNAME+"/", "constant/polyMesh/blockMeshDict") 
    parametrize_file(mesh_file, mesh_file, CFD_PARAMS)
    control_dict_file = os.path.join(CFD_FNAME+"/", "system/decomposeParDict") 
    parametrize_file(control_dict_file, control_dict_file, CFD_PARAMS)

    try:
        subprocess.check_output(["blockMesh", "-case", CFD_FNAME])
    except:
        assert False

    try:
        subprocess.check_output(["decomposePar", "-case", CFD_FNAME])
    except:
        assert False

    correct = run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
                       CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, True)
    if correct:
        # Reconstruct the fields from processor directories.
        try:
            subprocess.check_output(["reconstructPar", "-case", CFD_FNAME])
        except:
            assert False
        # Reconstruct the fields from processor directories.
        try:
            subprocess.check_output(["writeCellCentres", "-case", CFD_FNAME])
        except:
            assert False
        compare_forces(1e-6, CFD_PARAMS)
