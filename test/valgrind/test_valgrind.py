import pytest
from cplpy import run_test, prepare_config
import subprocess as sp
import os
from os.path import exists
from distutils.spawn import find_executable
import glob

import sys
udir = os.environ["CPL_PATH"] + '/test/'
sys.path.append(udir)
from testutils import cd, get_subprocess_error, runcmd

MD_EXEC = "./md"
CFD_EXEC = "./cfd"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))

def test_memory_leak(capsys):

    #Try to setup code
    buildstr = []
    with cd(TEST_DIR):
        buildstr.append(runcmd("rm -f md cfd"))
        buildstr.append(runcmd("cplf90 array_stuff.f90 md_recvsend_cells.f90 -o " + MD_EXEC))
        buildstr.append(runcmd("cplf90 array_stuff.f90 cfd_sendrecv_cells.f90 -o " + CFD_EXEC))

        with capsys.disabled():
            print("Builds info=", buildstr)

        #Check if executables have been built
        if (not exists(MD_EXEC)):
            raise IOError("Code md_recvsend_cells.f90 not compiling correctly")
        if (not exists(CFD_EXEC)):
            raise IOError("Code cfd_sendrecv_cells.f90 not compiling correctly")
        if (not find_executable("valgrind")):
            raise IOError("Valgrind not installed correctly")

        #Try to run code
        with capsys.disabled():
            print("running code")
    
        #Clean any previous runs
        runcmd("rm -f vg_*")

        #MPI based on repo mpich and valgrind raises a known error
        #cr_libinit.c:189 cri_init: sigaction() failed: Invalid argument
        #which can be safely ignored
        cmd = ("mpiexec -n 4 valgrind -v --leak-check=full --log-file='vg_md.%q{PMI_RANK}' " + MD_EXEC
                   + " : -n 2 valgrind -v --leak-check=full --log-file='vg_cfd.%q{PMI_RANK}' " + CFD_EXEC)
        out = runcmd(cmd)
        with capsys.disabled():
            print(out)

        #Check error
        filesgenerated = False
        files = glob.glob("vg_*")
        for filename in files:
            with open(filename,'r') as f:
                filestr = f.read()
                with capsys.disabled():
                    print("FILENAME=\n\n", filename, "\n\n\n\n", filestr)
                #Look to see if anything in definitely
                findstr= "definitely lost:"
                indx = filestr.find(findstr)
                if indx != -1:
                    line = filestr[indx+len(findstr):].split("\n")[0]
                    with capsys.disabled():
                        print("definitely lost = ", line)
                    assert int(line.split(" ")[1]) == 0
                else:
                    with capsys.disabled():
                        print("definitely lost not found, looking for other leak statements")
                    #If not, see if no leaks are possible statement
                    indx = filestr.find("no leaks are possible")
                    with capsys.disabled():
                        print("indx", indx)
                    assert indx != -1
                filesgenerated = True

    assert filesgenerated

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



