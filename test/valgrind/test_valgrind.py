import pytest
from cplpy import run_test, prepare_config
import subprocess as sp
import os
import glob

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_subprocess_error(e):
    print("subprocess ERROR")
    import json
    error = json.loads(e[7:])
    print((error['code'], error['message']))

# -----MAPPING TESTS-----

# EXPLANATION: These tests fail due to no_procs(MD) != k*no_procs(CFD),
#              k in [1,2,3,...] in one direction.

MD_EXEC = "./md"
CFD_EXEC = "./cfd"
TEST_TEMPLATE_DIR = os.path.join(os.environ["CPL_PATH"], "test/templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix():

    #Try to setup code
    mdcodes = "array_stuff.f90 md_recvsend_cells.f90"
    bldmd = ("mpif90 " + mdcodes 
             + "-I" + os.environ["CPL_PATH"] 
             + "/include  -L" + os.environ["CPL_PATH"] + "/lib  " 
             + "-Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./md")
    cfdcodes = "array_stuff.f90 cfd_sendrecv_cells.f90"
    bldcfd=  ("mpif90 " + cfdcodes 
             + " -I" + os.environ["CPL_PATH"] + "/include "
             + " -L" + os.environ["CPL_PATH"] + "/lib  " 
             + "-Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./cfd")
    with cd(TEST_DIR):
        try:
            out = sp.check_output("rm -f md cfd", shell=True)
            out = sp.check_output(bldmd, shell=True)
            out = sp.check_output(bldcfd, shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)


def test_memory_leak():

    #Try to run code
    cmd = ("mpiexec -n 4 valgrind --leak-check=full --log-file='vg_md.%q{PMI_RANK}' ./md "
               + ": -n 2 valgrind --leak-check=full --log-file='vg_cfd.%q{PMI_RANK}' ./cfd")
    with cd(TEST_DIR):
        try:
            out = sp.check_output("rm -f vg_*", shell=True)
            out = sp.check_output(cmd, shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith(b'error: {'):
                get_subprocess_error(e.output)


    #Check error
    files = glob.glob("vg_*")
    for filename in files:
        with open(filename,'r') as f:
            filestr = f.read()
            findstr= "definitely lost:"
            indx = filestr.find(findstr)
            line = filestr[indx+len(findstr):].split("\n")[0]
            print(line)
            assert int(line.split(" ")[1]) == 0

#@pytest.fixture()
#def prepare_config_fix(tmpdir):
#    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)
#    #Build code
#    try:
#        check_output("./build.sh", stderr=STDOUT, shell=True)
#    except:
#        raise

#@pytest.mark.parametrize("cfdprocs, mdprocs, err_msg", [
#                         ((1, 2, 1), (2, 2, 1), "")])
#def test_valgrind(prepare_config_fix, cfdprocs, mdprocs, err_msg):
#    MD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0}
#    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

#    CFD_PARAMS = {"lx": 24.0, "ly": 24.0, "lz": 24.0,
#                  "ncx": 24, "ncy": 24, "ncz": 24,
#                  "which_test": "cell_test"}
#    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

#    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
#                     "olap_xlo": 1, "olap_xhi": 24,
#                     "olap_ylo": 1, "olap_yhi": 4,
#                     "olap_zlo": 1, "olap_zhi": 24,
#                     "cnst_xlo": 1, "cnst_xhi": 1,
#                     "cnst_ylo": 1, "cnst_yhi": 1,
#                     "cnst_zlo": 1, "cnst_zhi": 1,
#                     "tstep_ratio": 50, }

#    parametrizeConfig(template_dir, config_params)



