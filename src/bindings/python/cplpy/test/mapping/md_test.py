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

if CPL.overlap():
    # Receiving cell coordinates from CFD
    olap_region = CPL.get_olap_limits()
    portion = CPL.my_proc_portion(olap_region)
    [ncxl, ncyl, nczl] = CPL.get_no_cells(portion)
    recv_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
    send_array = np.zeros((3, 0, 0, 0), order='F', dtype=np.float64)
    #CPL.scatter(send_array, olap_region, recv_array)
    CPL.recv(recv_array, olap_region)
    for imd in xrange(portion[0], portion[1] + 1):
        for jmd in xrange(portion[2], portion[3] + 1):
            for kmd in xrange(portion[4], portion[5] + 1):
                iloc, jloc, kloc, = CPL.map_glob2loc_cell(portion,
                                                          [imd, jmd, kmd])
                # Receive cell or coord depending on the test
                ixmd, jymd, kzmd = imd, jmd, kmd
                ixcfd, jycfd, kzcfd = recv_array[0:3, iloc, jloc, kloc]
                if which_test == "cell_test":
                    # This has to be true for every cell for the test to pass
                    ixcfd, jycfd, kzcfd = int(ixcfd), int(jycfd), int(kzcfd)
                elif which_test == "coord_test":
                    # This has to be true for every point for the test to pass
                    ixmd, jymd, kzmd = CPL.map_cell2coord(ixmd, jymd,
                                                          kzmd)
                    ixmd, jymd, kzmd = CPL.map_md2cfd_coord([ixmd, jymd,
                                                             kzmd])

                errix, errjy, errkz = ixcfd - ixmd, jycfd - jymd, kzcfd - kzmd
                if abs(errix) > 1e-6 or abs(errjy) > 1e-6 or abs(errkz) > 1e-6:
                    test_passed = False


                lines += str(ixcfd) + " " + str(jycfd) + " " + str(kzcfd) + " " +\
                         str(ixmd) + " " + str(jymd) + " " + str(kzmd) + "\n"


# If test do not pass save all the cell pairs to debug
lines = realm_comm.gather(lines, root=0)
myrank = realm_comm.Get_rank()
if myrank == 0:
    if which_test == "cell_test":
        fname = "md_recv_cells.dat"
    elif which_test == "coord_test":
        fname = "md_recv_coord.dat"
    with open(fname, "w") as cells_file:
        cells_file.writelines(lines)


if not test_passed:
    print("FAILED:", "There is something wrong in the mapping.",
          file=sys.stderr)
    comm_world.Abort(errorcode=1)


# Sending cell coordinates from MD to CFD
if CPL.overlap():
    send_array = np.zeros((3, ncxl, ncyl, nczl), order='F', dtype=np.float64)
    for iglob in xrange(portion[0], portion[1] + 1):
        for jglob in xrange(portion[2], portion[3] + 1):
            for kglob in xrange(portion[4], portion[5] + 1):
                iloc, jloc, kloc = CPL.map_glob2loc_cell(portion,
                                                         [iglob, jglob, kglob])
                if which_test == "cell_test":
                    send_array[0:3, iloc, jloc, kloc] = [iglob, jglob, kglob]
                elif which_test == "coord_test":
                    xglob, yglob, zglob = CPL.map_cell2coord(iglob, jglob,
                                                             kglob)
                    send_array[0:3, iloc, jloc, kloc] = [xglob, yglob, zglob]

    recv_array = np.zeros((3, 0, 0, 0), order='F', dtype=np.float64)
    CPL.send(send_array, olap_region)
