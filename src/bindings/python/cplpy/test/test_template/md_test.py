#!/usr/bin/env python2
from __future__ import print_function, division
from mpi4py import MPI
from cplpy import CPL
import numpy as np
import cPickle
import sys

comm_world = MPI.COMM_WORLD
CPL = CPL()

# Do not show any info to the stdin
CPL.set("output_mode", 0)

# Load parameters for the run
params = cPickle.load(open("md_params.dic", "rb"))

# Test selector flag
try:
    which_test = params["which_test"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)


# Parameters of the cpu topology (cartesian grid)
try:
    NPx = params["npx"]
    NPy = params["npy"]
    NPz = params["npz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)

NProcs = NPx*NPy*NPz
npxyz = np.array([NPx, NPy, NPz], order='F', dtype=np.int32)

# Parameters of the domain
try:
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)

xyzL = np.array([Lx, Ly, Lz], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.MD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processors is not coherent.", file=sys.stderr)
    comm_world.Abort(errorcode=1)

cart_comm = realm_comm.Create_cart([NPx, NPy, NPz])

CPL.setup_md(cart_comm, xyzL, xyz_orig)
lines = ""
test_passed = True

############################
#
# Write Tests Here
#
############################
if not test_passed:
    print("FAILED:", "Error message",
          file=sys.stderr)

