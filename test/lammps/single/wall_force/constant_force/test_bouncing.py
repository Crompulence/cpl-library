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
def clean_dir():

    print("Cleaning directory")
    #Try to setup code
    with cd(TEST_DIR):
        try:
            clean = sp.check_output("rm -f " + "./thermo_output* " 
                                             + "./log.lammps* " 
                                             + "./debug.vels" 
                                             + " " + MD_EXEC, shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)
            raise

    return clean

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
            raise

    return build

@pytest.fixture(scope="module")
def run_case():

    print("Running case")
    #Try to run code
    cmd = ('cplexec -m 1 "' + MD_EXEC + ' < single.in" ' + ' -c 1 ' +  CFD_EXEC)

    with cd(TEST_DIR):
        try:
            run = sp.check_output(cmd, shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)
            raise

    return run


@pytest.fixture(scope="module")
def build_run():
    try:
        clean = clean_dir()
        build = build_case()
    except sp.CalledProcessError:
        print("Build Failed")
    else:
        run = run_case()

def test_gravity(build_run):

    #Check vs analystical solution for gravity
    import bouncing

    with cd(TEST_DIR):
        error = bouncing.check_bouncing_error_vs_gravity()
        for e in error[0,1,:]:
            assert np.abs(e) < 1e-11


def test_regression(build_run):

    #Check vs analystical solution for gravity
    import bouncing as b
    with cd(TEST_DIR):

        t, z, v, f = b.read_data(logfile='./log.lammps', 
                                 datafile='./thermo_output.txt')

        t_reg, z_reg, v_reg, f_reg = b.read_data(logfile='./regression_data/log.lammps', 
                                     datafile='./regression_data/thermo_output.txt')

        for i in range(z.shape[0]):
            assert np.abs(z[i]-z_reg[i]) < 1e-12
            assert np.abs(v[i]-v_reg[i]) < 1e-12
            assert np.abs(f[i]-f_reg[i]) < 1e-12

