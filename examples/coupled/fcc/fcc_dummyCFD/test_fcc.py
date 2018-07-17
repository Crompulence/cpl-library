import numpy as np 

dp = 0.0565685
mu = 1.e-2
phi = 0.74
K = 435.
Ub = 0.01

eps = 1. - phi
Cd = 3*np.pi*mu*dp*K
Fy = Cd*Ub/eps

data = np.loadtxt('lammps/fcc200.dump', skiprows=9)
fy = data[:, 8]

print([np.min(fy), np.mean(fy), np.max(fy)])
print(Fy)

assert(np.all(abs(Fy - fy) <= 1e-6))