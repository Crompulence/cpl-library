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
from DragForce import DiFelice, Ergun

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
def get_input_parameters(md_input_file='./lammps/suzuki.in', cfd_case_dir='./openfoam/'):
    # LAMMPS input script
    mObj = LAMMPS_Input(md_input_file)

    # The porosity is hardwired to the value of a simple cubic array.
    mObj.epsf = 0.476401224402

    # OpenFOAM case 
    cObj = OpenFOAM_Input(case_dir=cfd_case_dir)
    _, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
    mObj.Uf = Uf[1]
    g = cObj.read_constant_file('./openfoam/constant/environmentalProperties', 'g')[1]
    rhof = cObj.read_constant_file('./openfoam/constant/transportProperties', 'rhob')
    nuf = cObj.read_constant_file('./openfoam/constant/transportProperties', 'nub')
    muf = nuf*rhof
    assert(mObj.gravity == -g)
    assert(mObj.fluid_density == rhof)
    assert(mObj.dynamic_viscosity == muf)

    return mObj

# Set the input parameters for the simulations. At present, only the inlet
# fluid velocity can be adjusted in the OpenFOAM input, while in the LAMMPS
# input, the drag model and movement type can be adjusted.
def set_input_parameters(Uf, dragModel, 
        md_input_file='./lammps/suzuki.in', cfd_case_dir='./openfoam/'):
    LAMMPS_Writer(md_input_file, 'dragModel', dragModel)
    OpenFOAM_Writer_0_File('./openfoam/0/Ub', 'inlet', Uf, isScalar=False, axis_val=1)

# Calculate the analytical pressure profile and the displacement of the top
# particle due to upward flow and buoyancy.
def analytical_pressure_displacement(mObj):
    
    dragModel = mObj.dragModel
    muf = mObj.dynamic_viscosity
    rhof = mObj.fluid_density
    Uf = mObj.Uf
    epsf = mObj.epsf
    dp = mObj.diameter
    rhop = mObj.density
    kn = mObj.kn
    ylo = mObj.ylo
    yhi = mObj.yhi
    s = mObj.lattice_scale
    g = -mObj.gravity 

    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero). Note that pSol is the
    # pressure at the outlet due to the drag force only and does not include
    # any contribution of the pressure gradient due to the presence of
    # gravity. This is added on when we get to the settlement calculation
    # below.
    if dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    else:
        raise ValueError('Unknown drag force model specified')

    gradP = fObj.calculate_pressure_gradient(Uf=Uf/epsf, Vp=0.)
    H = yhi-ylo
    pSol = gradP*H

    fdi = fObj.calculate_drag_force(Uf=Uf/epsf, Vp=0.)
    volume = (np.pi/6)*dp**3
    fi = volume*rhop*g - volume*rhof*g + fdi
    N = H/s
    disp = 0.5*N*(N+1)*fi/kn
    xySol = (yhi - 0.5*dp) + disp
    
    return pSol, xySol

# Read print data for the top particle on column
def read_print_data(xy0, print_file='./lammps/print_suzuki.txt'):
    data = np.loadtxt(print_file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    # Append initial values
    t = np.insert(t, 0, 0)
    xy = np.insert(xy, 0, xy0)

    return t, xy

# Obtain the pressure profile from the numerical simulation
def get_pressure_profile(CFD_folder='./openfoam/', pressure_var='p', axis_val=1):
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    pObj = PPObj.plotlist[pressure_var]
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    p = np.squeeze(p)

    return h, p

# Plot pressure profile and displacement obtained from numerical simulation
# and analytical solution. Save the file in the results directory (which is
# created if required) and also save the data used for plotting as .npz file.
def plot_pressure_displacement(t, xy, xySol, h, p, pSol, ylo, yhi, file_name='./fig_pressure'):
    # Plot pressure
    plt.plot(h, p, 'r-')
    plt.plot([ylo,yhi], [pSol,0.], 'k--')
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
    plt.savefig(file_name + '_pressure.png')
    np.savez(file_name + '_pressure.npz', h=h, p=p, pSol=pSol, ylo=ylo, yhi=yhi)
    plt.close()

    # Plot displacement
    plt.plot(t, xy, 'r-')
    plt.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.xlabel('Time (s)')
    plt.ylabel('Position of Top Particle (cm)')
    plt.legend(('Numerical', 'Analytical (Upward Flow)'))
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '_displacement.png')
    np.savez(file_name + '_displacement.npz', t=t, xy=xy, xySol=xySol)
    plt.close()
    

# Compare the final displacement of the top-most particle at the end of the
# simulation for a specified relative error.
def compare_displacement(xy, xySol, tol=0.01):
    err = abs((xySol -xy[-1])/xySol) <= tol
    assert err, ('Final displacement of {:.6f} does not match analytical'.format(xy[-1])
                 +' solution of {:.6f} within {:.2f}% relative error.'.format(xySol, tol*100))

# Compare the pressure along the length of the sample at the end of the
# simulation for a specified relative error.
def compare_pressure(h, p, pSol, ylo, yhi, tol=0.01):
    pSolProfile = np.interp(h, [ylo,yhi], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = abs((pSolProfile[i] - p[i])/pSolProfile[i]) <= tol
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     +' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))

# ----- Main ----- #
dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
@pytest.mark.parametrize('Uf', Uf_values)
@pytest.mark.parametrize('dragModel', dragModels)
def test_displacement(Uf, dragModel):

    # Set input parameters
    set_input_parameters(Uf, dragModel)

    # Run coupled simulation
    run_coupled()
    if not os.path.exists('lammps/print_suzuki.txt'):
        print('Attempting to re-run coupled run.')
        run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Load print data
    xy0 = mObj.yhi - 0.5*mObj.diameter
    t, xy = read_print_data(xy0)

    # Get pressure profile
    h, p = get_pressure_profile()
    p = p - mObj.fluid_density*mObj.gravity*np.flipud(h)

    # Extract input parameters from lammps input script
    pSol, xySol = analytical_pressure_displacement(mObj)

    # Plot the results
    plot_pressure_displacement(t, xy, xySol, h, p, pSol, mObj.ylo, mObj.yhi,
        file_name='./results/fig_Uf_{}_{}'.format(Uf, dragModel))

    # Test analytical solution
    compare_displacement(xy, xySol, tol=0.02)
    if Uf > 0.0:
        compare_pressure(h, p, pSol, mObj.ylo, mObj.yhi, tol=0.02)
    
if __name__ == "__main__":
    test_displacement(Uf=0.0, dragModel='DiFelice')