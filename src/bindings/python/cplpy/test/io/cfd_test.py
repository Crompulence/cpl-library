#!/usr/bin/env python2
from __future__ import print_function, division
from mpi4py import MPI
from cplpy import CPL, CPL_VAR_TYPES
import numpy as np
import cPickle
import sys
import json

comm_world = MPI.COMM_WORLD
CPL = CPL()

# Do not show any info to the stdin
CPL.set("output_mode", 0)

# Load parameters for the run
params = cPickle.load(open("cfd_params.dic", "rb"))

# Parameters of the cpu topology (cartesian grid)
try:
    NPx = params["npx"]
    NPy = params["npy"]
    NPz = params["npz"]
except:
    print("ERROR: ", sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)


# NPx = 3
# NPy = 1
# NPz = 3

NProcs = NPx*NPy*NPz

# Create communicators and check that number of processors is consistent
realm_comm = CPL.init(CPL.CFD_REALM)
nprocs_realm = realm_comm.Get_size()

if (nprocs_realm != NProcs):
    print("ERROR: ", "Number of processors is not coherent.", file=sys.stderr)
    comm_world.Abort(errorcode=1)

json_fname = "config.cpl"
json_nocomments_fname = "config_nocomments.cpl"
file_dict = {}
local_dict = json.load(open(json_nocomments_fname))
CPL.load_param_file(json_fname)
TYPES = CPL_VAR_TYPES
fields = {"string": TYPES.STRING, "int": TYPES.INT, "int_array": TYPES.INT_ARRAY,
          "double": TYPES.DOUBLE, "double_array": TYPES.DOUBLE_ARRAY, "bool_true": TYPES.BOOL,
          "bool_false": TYPES.BOOL}
sections = ["md", "cfd"]
subsections = ["", "subsection"]
myrank = realm_comm.Get_rank()
for s in sections:
    file_dict[s] = {}
    for sub in subsections:
        if sub:
            file_dict[s][sub] = {}
        for f, ftype in fields.items():
            if sub:
                file_dict[s][sub][f] = CPL.get_file_var(s, f, ftype)
            else:
                file_dict[s][f] = CPL.get_file_var(s, f, ftype)
if local_dict == file_dict:
    print ("CFD proc %d: SUCCESS!" % myrank)
else:
    print("CFD proc %d: FAILURE!" % myrank, sys.exc_info()[0], file=sys.stderr)
    comm_world.Abort(errorcode=1)
CPL.finalize()
