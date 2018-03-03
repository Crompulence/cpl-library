import numpy as np
import subprocess as sp
import os
import pytest

plotstuff=False

#Compare a range of cases
rho = 1e3
mu = 0.001

outputdir="./drag_output"

#File and running stuff
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

#Generate drag results from CPL_force
print("Rebuilding and running c++ code")
with cd("../../"):

    try:
        cplbuild = sp.check_output("make", shell=True)
    except sp.CalledProcessError as e:
        if e.output.startswith('error: {'):
            get_subprocess_error(e.output)
        raise

try:
    dragtestbuild = sp.check_output("make", shell=True)
    mkdir_output = sp.check_output("mkdir -p " + outputdir, shell=True)
    dragtestrun = sp.check_output("./run_CPLForce_drag_models", shell=True)
except sp.CalledProcessError as e:
    if e.output.startswith('error: {'):
        get_subprocess_error(e.output)
    raise

import sys

#Get drag utils from Chris' work
dragutilsdir = "./drag-utils"
download = False
if os.path.exists(dragutilsdir):
    try:
        sys.path.insert(0, dragutilsdir+"/")
        import drag_utils.correlations as c
    except ImportError:
        download=True
        os.rename(dragutilsdir, dragutilsdir+".bak")
else:
    download = True

if download:
    #Change this back to Chris's branch once he excepts pull request
    print("Downloading Chris Knight's drag-utils")
    getdragutils = sp.check_output("git clone https://github.com/edwardsmith999/drag-utils.git "+dragutilsdir, shell=True)
    sys.path.insert(0, dragutilsdir)
    import drag_utils.correlations as c


def get_data(case):

    #Load CPLForce data from file
    data = np.genfromtxt(outputdir+"/"+case, delimiter=",", names=True)

    #Get force type from Chris' library
    Forcefn = getattr(c, case.replace("_",""))
    Fpy = Forcefn(data['phi'], data['D'], data['v0'], rho=rho, mu=mu, norm=False)

    return data, Fpy

@pytest.mark.parametrize("case,out", [
    ("Di_Felice", True),
    ("Stokes", True),
    ("Ergun", True),
    ("BVK", True),
    ("Tang", True),
])
def test_answer(case,out):
    data, Fpy = get_data(case)
    print(case + " max Error = ", np.max(np.abs((Fpy+data['F0'])/Fpy)))
    assert(np.max(np.abs((Fpy+data['F0'])/Fpy) < 1e-5) == out)


cases = ["Di_Felice", "Stokes", "Ergun", "BVK", "Tang"]

if plotstuff:
    import matplotlib.pyplot as plt
    from matplotlib import rc

    #Plotting stuff
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    for case in cases:

        CPLdata, Fpy = get_data(case)

        #Plot Both
        plt.plot(CPLdata['phi'], -Fpy, 'k-', label="Chris' Python Script")
        plt.plot(CPLdata['phi'], CPLdata['F0'], 'ro', label="CPLForce")
        #plt.plot(CPLdata['phi'], (Fpy+CPLdata['F0'])/Fpy, 'b-', label="Error")

        plt.legend(loc=3)
        plt.title(case)
        plt.xlabel(r"$\phi$")
        plt.ylabel("$F$")
        plt.show()
        #plt.savefig(case + "phi.pdf", bbox_inches="tight")

        plt.plot(CPLdata['D'], -Fpy, 'k-', label="Chris' Python Script")
        plt.plot(CPLdata['D'], CPLdata['F0'], 'ro', label="CPLForce")
        #plt.plot(CPLdata['D'], (Fpy+CPLdata['F0'])/Fpy, 'b-', label="Error")

        plt.legend(loc=3)
        plt.title(case)
        plt.xlabel(r"$D$")
        plt.ylabel("$F$")
        plt.show()

        print(case + " max Error = ", np.max(np.abs((Fpy+CPLdata['F0'])/Fpy)))
        assert(np.max(np.abs((Fpy+CPLdata['F0'])/Fpy) < 1e-5) == True)
        #plt.savefig(case + "D.pdf", bbox_inches="tight")

    
