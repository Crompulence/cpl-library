import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

# Add python scripts to path and import required classes
sys.path.append('../../python_scripts/')
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

# Extract the input parameters from DEM script for LAMMPS and OpenFOAM case
# directory. Check that they are consistent.
def get_input_parameters(md_input_file='./lammps/resting.in', cfd_case_dir='./openfoam'):    
    # LAMMPS script
    mObj = LAMMPS_Input(md_input_file)

    # Note cell volume hard-wired (for ease)
    CELL_VOLUME = 0.1**3
    mObj.epsf = (CELL_VOLUME - (np.pi/6)*mObj.diameter**3)/CELL_VOLUME

    cObj = OpenFOAM_Input(case_dir=cfd_case_dir)
    _, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
    mObj.Uf = Uf[1]

    return mObj

# Set the input parameters for the simulations. At present, only the particle
# diameter and drag force model can be adjusted. Both these only apply to the
# LAMMPS input.
def set_input_parameters(dp, y0, dragModel, Uf, 
        md_input_file='./lammps/resting.in', cfd_case_dir='./openfoam/'):
    
    LAMMPS_Writer(md_input_file, 'diameter', dp)
    LAMMPS_Writer(md_input_file, 'y0', y0)
    LAMMPS_Writer(md_input_file, 'dragModel', dragModel)
    OpenFOAM_Writer_0_File('./openfoam/0/Ub', 'inlet', Uf, isScalar=False, axis_val=1)

# Calculate the analytical displacement and velocity profile of the particle,
# along with the terminal velocity. This is only applicable for the Stokes
# drag law, at present.
def analytical_displacement(mObj):

    dragModel = mObj.dragModel
    muf = mObj.dynamic_viscosity
    rhof = mObj.fluid_density
    Uf = mObj.Uf
    epsf = mObj.epsf
    dp = mObj.diameter
    rhop = mObj.density
    kn = mObj.kn
    y0 = mObj.y0
    g = -mObj.gravity 

    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    else:
        raise('Unknown drag force model specified')

    if y0 < 1.0:
        CORRECTION_FACTOR = 0.5
    else:
        CORRECTION_FACTOR = 1.0

    fdi = fObj.calculate_drag_force(Uf=Uf/epsf, Vp=0.)
    print(fdi)
    volume = (np.pi/6)*dp**3
    fi = volume*rhop*g - CORRECTION_FACTOR*volume*rhof*g + CORRECTION_FACTOR*fdi
    disp = fi/kn
    xySol = y0 + disp
    
    return xySol

# Read print data for the top particle on column
def read_print_data(xy0, print_file='./lammps/print_resting.txt'):
    data = np.loadtxt(print_file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    # Append initial values
    t = np.insert(t, 0, 0)
    xy = np.insert(xy, 0, xy0)

    return t, xy

# Plot displacement and velocity profile obtained from numerical simulation
# and analytical solution. Save the file in the results directory (which is
# created if required) and also save the data used for plotting as .npz file.
def plot_displacement(t, xy, xySol, file_name='./fig'):
    # Plot displacement
    plt.plot(t, xy, 'r-')
    plt.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.xlabel('Time (s)')
    plt.ylabel('Position (cm)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    plt.ticklabel_format(useOffset=False)
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '.png')
    np.savez(file_name + '.npz', t=t, xy=xy, xySol=xySol)
    plt.close()

# Compare the displacement profile with time for a specified relative error.
def compare_displacement(xy, xySol, tol=0.01):
    err = abs((xySol - xy[-1])/xySol) <= tol
    assert err, ('Displacement of {:.6f} does not match analytical'.format(xy[-1])
            + ' solution of {:.6f} within {:.2f}% relative error.'.format(xySol, tol*100))

# ----- Main ----- #
# dragModels = ['Drag', 'Stokes', 'DiFelice', 'Ergun']
dragModels = ['DiFelice']
dp_values = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
y0_values = [1.0]
Uf_values = [40.]
@pytest.mark.parametrize('Uf', Uf_values)
@pytest.mark.parametrize('dragModel', dragModels)
@pytest.mark.parametrize('y0', y0_values)
@pytest.mark.parametrize('dp', dp_values)
def test_displacement(dp, y0, dragModel, Uf):

    # Set input parameters
    y0ini = y0
    y0 = y0 + 0.5*dp
    set_input_parameters(dp, y0, dragModel, Uf)

    # Run coupled simulation
    run_coupled()
    if not os.path.exists('lammps/print_resting.txt'):
        print('Attempting to re-run coupled run.')
        run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Load print data
    t, xy = read_print_data(mObj.y0)
    
    # Extract input parameters from lammps input script
    xySol = analytical_displacement(mObj)

    # Plot the results
    plot_displacement(t, xy, xySol,
        file_name='./results/fig_{}_dp_{}_Uf_{}_y0_{:.1f}'.format(dragModel, dp, Uf, y0ini))

    # Test analytical solution
    compare_displacement(xy, xySol, tol=0.01)
    
if __name__ == "__main__":
    test_displacement(dp=0.10, y0=1.00, dragModel='DiFelice', Uf=40.0)
