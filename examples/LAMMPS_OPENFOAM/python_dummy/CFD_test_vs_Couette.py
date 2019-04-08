############################################################################
#
#    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
#     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
#      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
#       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
#        _\/\\\_____________\/\\\/////////____\/\\\_____________
#         _\//\\\____________\/\\\_____________\/\\\_____________
#          __\///\\\__________\/\\\_____________\/\\\_____________
#           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
#            _______\/////////__\///______________\///////////////__
#
#
#                         C P L  -  L I B R A R Y
#
#           Copyright (C) 2012-2018 Edward Smith & David Trevelyan
#
# License
#
#    This file is part of CPL-Library.
#
#    CPL-Library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CPL-Library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.
#
# Description
#
#    Minimal CFD MOCK Script for debugging
#
# Author(s)
#
#    Edward Smith
#


import numpy as np
from mpi4py import MPI

from CouetteAnalytical import CouetteAnalytical as CA

from cplpy import CPL

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)

## Parameters of the cpu topology (cartesian grid)
npxyz = [1, 1, 1]

# FCC lattice, plus wall width plus buffer region
xyzL = [16.795961913825074, 45.349097, 16.795961913825074]
xyz_orig = [0.0, 0.0, 0.0]
ncxyz = [8, 8, 8]

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=3)

#Get values from CPL library
ncy = CPL.get("ncy")
yL_md = CPL.get("yl_md")
yL_cfd = CPL.get("yl_cfd")
dy = CPL.get("dy")
jcmax_olap = CPL.get("jcmax_olap")

#Setup Coupled domain details
yL_cpl = yL_md + yL_cfd - dy*jcmax_olap
ncy_cpl = int(yL_cpl/dy)

#Create analytical solution object
Uwall = 1.0
nu = 1.7
Re = yL_cpl/nu
CAObj = CA(Re=Re, U=Uwall, Lmin=0.0, Lmax=yL_cpl, npoints=ncy_cpl+2)


for time in range(500):

    y_anal, u_anal = CAObj.get_vprofile(time)

    # send data to update
    send_array[...] = 0.0
    send_array[0,:,:,:] = Uwall # u_anal[-ncy]
    CPL.send(send_array)
        
    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array)
    Ux = np.mean(recv_array[0,:,:,:])/np.mean(recv_array[3,:,:,:])

    print("CFD time = ", time, Ux, u_anal[-ncy])
#    import matplotlib.pyplot as plt
#    plt.plot(y_anal, u_anal)
#    plt.plot(y_anal[-ncy], u_anal[-ncy], 'ro')
#    plt.show()

comm.Barrier()
CPL.finalize()
MPI.Finalize()




