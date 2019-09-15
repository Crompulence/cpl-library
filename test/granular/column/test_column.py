import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

# Import postproclib
sys.path.insert(0, "./pyDataView/")
try:
    import postproclib as ppl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input, LAMMPS_Writer
from OpenFOAM_Input import OpenFOAM_Input, OpenFOAM_Writer_0_File
from DragForce import DragForce, Stokes, DiFelice, Ergun

# Run coupled simulation as subprocess
def run_coupled(run_bash_script='run.sh'):
    try:
        cmd = './' + run_bash_script
        p = sp.Popen(cmd, 
            stdout=sp.PIPE, 
            stderr=sp.STDOUT, 
            shell=True)
        while p.poll() is None:
            l = p.stdout.readline()
            print(l.rstrip())
        print(p.stdout.read())
    except:
        raise RuntimeError('Error running bash run script' + run_bash_script + ' in base directory.')

# Read print data for the top particle on column
def read_print_data(print_file):
    data = np.loadtxt(print_file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    return t, xy

# Extract the input parameters from DEM script for LAMMPS and OpenFOAM case
# directory. Check that they are consistent.
def get_input_parameters(md_input_file='./lammps/column.in', cfd_case_dir='./openfoam/'):    
    # LAMMPS script
    mObj = LAMMPS_Input(md_input_file)
    mObj.xyz_orig = np.array([mObj.xlo, mObj.ylo, mObj.zlo], order='F', dtype=np.float64)
    mObj.xyzL = np.array([mObj.xhi, mObj.yhi, mObj.zhi], order='F', dtype=np.float64)

    # The porosity is hardwired to the value of a simple cubic array.
    mObj.epsf = 0.476401224402

    # OpenFOAM case 
    cObj = OpenFOAM_Input(cfd_case_dir)
    _, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
    mObj.Uf = Uf[1]
    rhof = cObj.read_constant_file('./openfoam/constant/transportProperties', 'rhob')
    nuf = cObj.read_constant_file('./openfoam/constant/transportProperties', 'nub')
    muf = nuf*rhof
    assert(mObj.fluid_density == rhof)
    assert(mObj.dynamic_viscosity == muf)

    return mObj

# Set the input parameters for the simulations. At present, only the inlet
# fluid velocity can be adjusted in the OpenFOAM input, while in the LAMMPS
# input, the drag model and movement type can be adjusted.
def set_input_parameters(Uf, dragModel, movementType, 
        md_input_file='./lammps/column.in', cfd_case_dir='./openfoam/'):
    
    LAMMPS_Writer(md_input_file, 'dragModel', dragModel)
    LAMMPS_Writer(md_input_file, 'movementType', movementType)
    OpenFOAM_Writer_0_File('./openfoam/0/Ub', 'inlet', Uf, isScalar=False, axis_val=1)

# Calculate the analytical solution for the pressure at the inlet, assuming
# that the pressure at the outlet is zero
def analytical_pressure(mObj):

    dragModel = mObj.dragModel
    muf = mObj.dynamic_viscosity
    rhof = mObj.fluid_density
    Uf = mObj.Uf
    dp = mObj.diameter
    epsf = mObj.epsf
    xyzL = mObj.xyzL
    xyz_orig = mObj.xyz_orig

    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero)
    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    else:
        raise('Unknown drag force model specified')

    gradP = fObj.calculate_pressure_gradient(Uf=Uf/epsf, Vp=0.)
    pSol = gradP*(xyzL[1] - xyz_orig[1])

    return pSol
    
# Obtain the pressure profile from the numerical simulation
def get_pressure_profile(CFD_folder='./openfoam/', pressure_var='p', axis_val=1):
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    pObj = PPObj.plotlist[pressure_var]
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    p = np.squeeze(p)

    return h, p

# Plot pressure profile obtained from numerical simulation and analytical
# solution. Save the file in the results directory (which is created if
# required) and also save the data used for plotting as .npz file.
def plot_pressure(h, p, pSol, xyz_orig, xyzL, file_name='fig_pressure'):
    plt.plot(h, p, 'r-')
    plt.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (0.1Pa)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '.png')
    np.savez(file_name + '.npz', h=h, p=p, pSol=pSol, xyz_orig=xyz_orig, xyzL=xyzL)
    plt.close()

# Compare the pressure along the length of the sample at the end of the
# simulation for a specified relative error.
def compare_pressure(h, p, pSol, xyz_orig, xyzL, tol=0.01):
    pSolProfile = np.interp(h, [xyz_orig[1],xyzL[1]], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = abs((pSolProfile[i] - p[i])/pSolProfile[i] <= tol)
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     +' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))

# ----- Main ----- #
movementTypes = ['fixed', 'moving']
dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.1, 0.2, 0.3, 0.4, 0.5]
@pytest.mark.parametrize('movementType', movementTypes)
@pytest.mark.parametrize('dragModel', dragModels)
@pytest.mark.parametrize('Uf', Uf_values)
def test_pressure(Uf, dragModel, movementType):
    # Set input parametesr
    set_input_parameters(Uf, dragModel, movementType)

    # Run coupled simulation
    run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Calculate analytical pressure profile
    pSol = analytical_pressure(mObj)

    # Extract numerical pressure profile (end of simulation)
    h, p = get_pressure_profile()

    # Plot the results
    plot_pressure(h, p, pSol, mObj.xyz_orig, mObj.xyzL, 
        file_name='./results/fig_pressure_Uf_{}_{}_{}'.format(Uf, dragModel, movementType))

    compare_pressure(h, p, pSol, mObj.xyz_orig, mObj.xyzL, tol=0.02)

if __name__ == "__main__":
    test_pressure(Uf=0.5, dragModel='DiFelice', movementType='fixed')
