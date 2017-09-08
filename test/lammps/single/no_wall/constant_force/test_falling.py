import pytest
from cplpy import run_test, prepare_config
import subprocess as sp
import os
import glob
import numpy as np

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
    print(error['code'], error['message'])


MD_EXEC = "./lmp_cpl"
CFD_EXEC = "./CFD_single_ball.py"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))

@pytest.fixture(scope="module")
def build_case():
    print("Building LAMMPS")
    #Try to setup code
    with cd(TEST_DIR):
        try:
            build = sp.check_output("./build.sh", shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)

    return build

@pytest.fixture(scope="module",
                params=[2., 5., 9.81, 14.])
def run_case(request):

    #Try to run code
    cmd = ('cplexec -m 1 "' + MD_EXEC + ' < single.in" ' + ' -c 1 "' +  CFD_EXEC + " " + str(request.param) + ' "')

    print("Running case " + cmd)
    with cd(TEST_DIR):
        try:
            clean = sp.check_output("rm -f ./thermo_output* ./log.lammps* ./debug.vels", 
                                    shell=True)
            run = sp.check_output(cmd, shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)

    return request


#@pytest.fixture(scope="module")
#def build_run():
#    build = build_case()
#    run = run_case()

def test_gravity(build_case, run_case):

    #Check vs analystical solution for gravity
    import falling

    with cd(TEST_DIR):
        error = falling.check_falling_error_vs_gravity(g=run_case.param)
        for e in error:
            assert np.abs(e) < 1e-14


