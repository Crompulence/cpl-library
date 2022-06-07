import pytest
import cplpy
from cplpy import run_test, prepare_config
import subprocess as sp
import os
import glob
import time

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


def runcmd(cmd):
    try:
        p = sp.Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        output, error = p.communicate()
        if p.returncode != 0: 
            print("returncode", p.returncode)
            print("Stdout = ", output)
            print("Stderror = ", error)
        #run = sp.check_output(cmd, stderr=sp.STDOUT, shell=True).decode("utf-8")
    except sp.CalledProcessError as e:
        print("Stdout = ", e.stdout)
        print("Stderror = ", e.stderr)
        if e.output.startswith(b'error: {'):
            get_subprocess_error(e.output)
        raise

    return run

if (cplpy.CPL.MPI_version == "OPENMPI"):
    cases=["split"]
else:
    cases=["port", "split"]
perms = []
dirs = ["fortran", "cpp", "python"]
for cplt in cases:
    for mdd in dirs:
        for cfdd in dirs:
            perms.append([cplt, mdd, cfdd, 8])

@pytest.mark.parametrize("cpltype, mddir, cfddir, count", perms)
def test_example_sendrecv_scripts(cpltype, mddir, cfddir, count):

    #Try to run code
    CDIR = os.path.dirname(os.path.realpath(__file__))
    EXAMPLES_DIR =  os.path.join(CDIR, "../../examples/sendrecv_globcell/")

    with cd(EXAMPLES_DIR):

        if (cpltype is "split"):
            if (cplpy.CPL.MPI_version == "OPENMPI"):
                cmd = "mpiexec --oversubscribe -n 16 "
            else:
                cmd = "mpiexec -n 16 "
        elif cpltype is "port":
            cmd = "cplexec -m 16 "

        if ("python" in mddir):
            cmd += "python ./" + mddir + "/md_recv_cells.py"
        else:
            with cd(mddir):
                runcmd("./build.sh")
            cmd += "./" + mddir + "/md"

        #Choose between split 
        if cpltype is "split":
            cmd += " : -n 4 "
        elif cpltype is "port":
            cmd += " -c 4 "
        else:
            raise ValueError("cpltype type unknown", cpltype)

        if ("python" in cfddir):
            cmd += " python ./" + cfddir + "/cfd_send_cells.py"
        else:
            with cd(cfddir):
                runcmd("./build.sh")
            cmd += " ./" + cfddir + "/cfd"

        if cpltype is "port":
            cmd = cmd.replace("python ", "")

        print(cmd)
        runstr = runcmd(cmd)

        #Check for 8 recv counts
        c = 0
        findstr= "CELLS HAVE BEEN RECEIVED CORRECTLY"
        for l in runstr.split("\n"):
            if (l.find(findstr) != -1):
                c += 1

        assert c == count

