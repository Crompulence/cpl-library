import os
import sys
import subprocess as sp
import numpy as np
import pytest

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
from LAMMPS_Input import LAMMPS_Input
from DragForce import DragForce, Stokes, DiFelice

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
        print('Error running bash run script' + run_bash_script + ' in base directory.')
        raise

def read_print_data(file):
    data = np.loadtxt(file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    return t, xy

def analytical_solution_from_input_parameter(file):
    mdObj = LAMMPS_Input();
    muf = mdObj.read_variables(file, 'dynamic_viscosity')
    rhof = mdObj.read_variables(file, 'fluid_density')
    Uf = mdObj.read_variables(file, 'fluid_velocity')
    dp = mdObj.read_variables(file, 'diameter')
    rhop = mdObj.read_variables(file, 'density')
    kn = mdObj.read_variables(file, 'kn')
    g = mdObj.read_variables(file, 'gravity')
    s = mdObj.read_variables(file, 'lattice_scale')
    yhi = mdObj.read_variables(file, 'yhi')
    ylo = mdObj.read_variables(file, 'ylo')
    dragModel = mdObj.read_variables(file, 'dragModel')
    
    # The porosity is hard-wired because the configuration for the Suzuki
    # column test is a simple cubic arrangement.
    epsf = 0.4764

    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero)
    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=0., epsf=epsf)
    else:
        raise('Unknown drag force model specified')

    # Add pressure gradient contribution from gravity and calculate pressure
    # at the inlet
    gradP = fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf, Vp=0., dp=dp)
    H0 = yhi-ylo
    pSol = (gradP + rhof*g)*H0

    # Obtain the swelling displacement due to upward seepage flow. Initially,
    # calculate the coefficient of volume compressibility assuming a
    # settlement due to gravity and buoyancy (even though gravity or buoyancy
    # is not simulated).
    rhosat = (1-epsf)*rhop + epsf*rhof
    sigma = (rhosat-rhof)*g*(H0/2)
    fi = 4./3.*np.pi*(dp/2)**3*(rhop-rhof)*g
    N = H0/s
    settlement = 0.5*N*(N+1)*fi/kn
    mv = settlement/(sigma*H0)

    # Height after settlement due to gravity and buoyancy. Also the position
    # of top most particle (assumed to be radius away from top boundary) after
    # settlement.
    H = H0 - settlement
    xySol_settle = yhi - 0.5*dp - settlement

    # Swelling due to upward flow and position of top most particle after
    # upward flow. 
    disp = 0.5*mv*(H**2)*gradP
    xySol_upward = xySol_settle + disp
    
    return Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward
    
def get_velocity_pressure_porosity_profile(CFD_folder='./openfoam/', 
    velocity_var='Ub', pressure_var='p', porosity_var='alpha', axis_val=1):
    # Get plotting object
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    UObj = PPObj.plotlist[velocity_var]
    pObj = PPObj.plotlist[pressure_var]
    nObj = PPObj.plotlist[porosity_var]

    # Get profile (note that CPL calculates 1-eps)
    h, U = UObj.profile(axis=axis_val, startrec=UObj.maxrec, endrec=UObj.maxrec)
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    h, n = nObj.profile(axis=axis_val, startrec=nObj.maxrec, endrec=nObj.maxrec)
    h = h
    U = U[:,1]
    p = np.squeeze(p)
    n = np.squeeze(np.ones_like(n)-n)

    return h, U, p, n

def compare_displacement(xySol_upward, xy, tol):
    # Test final displacement matches analytical solution.
    err = abs((xySol_upward -xy[-1])/xySol_upward <= tol)
    assert err, ('Final displacement of {:.6f} does not match analytical'.format(xy[-1])
                 +' solution of {:.6f} within {:.2f}% relative error.'.format(xySol_upward, tol*100))

def compare_pressure(pSol, p, h, ylo, yhi, tol=0.01):
    # Test pressure profile at the end of simulation.
    pSolProfile = np.interp(h, [ylo,yhi], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = abs((pSolProfile[i] - p[i])/pSolProfile[i] <= tol)
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     +' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))


@pytest.fixture(scope="module")
def setup():

    # Run coupled simulation
    run_coupled()

    # Load print data
    t, xy = read_print_data('lammps/print_column.txt')

    # Extract input parameters from lammps input script
    Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward = analytical_solution_from_input_parameter('lammps/column.in')

    # Extract velocity, pressure and eps profile (at end of simulation only)
    h, U, p, n = get_velocity_pressure_porosity_profile()

    return t, xy, Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward, h, U, p, n


tols = [1.0, 0.1, 0.01, 0.001, 0.0002]
@pytest.mark.parametrize("tols", tols)
def test_pressure(setup, tols):
    t, xy, Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward, h, U, p, n = setup
    compare_pressure(pSol, p, h, ylo, yhi, tols)

tols = [1.0, 0.1, 0.01, 0.001, 0.0001, 9e-5]
@pytest.mark.parametrize("tols", tols)
def test_displacement(setup, tols):
    t, xy, Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward, h, U, p, n = setup
    compare_displacement(xySol_upward, xy, tols)


if __name__ == "__main__":

    t, xy, Uf, epsf, pSol, ylo, yhi, xySol_settle, xySol_upward, h, U, p, n = setup()

    import matplotlib.pyplot as plt

    # Plot Displacement Profile
    plt.plot(t, xy, 'r-')
    plt.plot(t, np.ones_like(t)*xySol_upward, 'k--')
    # plt.plot(t, np.ones_like(t)*xySol_settle, 'k:')
    plt.xlabel('Time (s)')
    plt.ylabel('Position of Top Particle (cm)')
    plt.legend(('Numerical', 'Analytical (Upward Flow)'))
    plt.tight_layout()
    plt.savefig('fig_displacement.png')
    plt.close()

    # Plot Pressure Profile
    plt.plot(h, p, 'r-')
    plt.plot([ylo,yhi], [pSol,0.], 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (0.1Pa)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    plt.savefig('fig_pressure.png')
    plt.close()
