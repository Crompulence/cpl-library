import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def analytical_gravity(z0, v0, t0, t, g=9.81):
    return z0 + v0*(t-t0) - 0.5*g*(t-t0)**2


def read_data(logfile='./log.lammps', 
              datafile='./thermo_output.txt'):

    #Get data from file
    with open(logfile) as f:
        for l in f.readlines():
            if l.find("timestep") != -1:
                dt = float(l.strip('timestep'))
                break
        
    #print("timestep = ", dt)
    data = np.genfromtxt(datafile)
    t = data[:,0]*dt
    z = data[:,1]
    v = data[:,2]
    f = data[:,3]

    return t, z, v, f

def unwrap(z):

    Lz = 2.5e-3    
    znew = np.zeros(z.shape[0])
    znew[0] = z[0]
    for i in range(1,z.shape[0]):
        dz = z[i] - z[i-1]
        if dz > 0.:
            dz -= Lz
        znew[i] = znew[i-1] + dz
    
    return znew


def check_falling_error_vs_gravity(D=3.5e-4, g=9.81):

    t, z, v, f = read_data()
    
    #Error vs gravity
    za = analytical_gravity(z[0], v[0], t[0], t, g)
    error = unwrap(z) - za

    return np.array(error)


if __name__ == "__main__":
    
    error = check_falling_error_vs_gravity()

    plt.plot(error)
    plt.show()
