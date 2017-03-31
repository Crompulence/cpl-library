import numpy as np
import sys


def check_forces(tol, cfd_fname="cfd_forces.dat", lammps_fname="lammps_forces.dat"):

    # Line format of CFD script file -- > x z vx vy vz
    with open(cfd_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
    cfd_cells = {}
    for l in cfd_lines:
        cfd_cells[(float(l[0]), float(l[1]), float(l[2]))] = np.array([float(l[3]),
                                                                       float(l[4]),
                                                                       float(l[5])])

    # Line format of LAMMPS file -- > chunk x y z ncount fx fy fz
    with open(lammps_fname, "r") as lammps_file:
        lammps_lines = lammps_file.readlines()
    skip = int(lammps_lines[3].split(" ")[1])
    lammps_lines = lammps_lines[4:]
    lammps_lines = lammps_lines[skip+1:]
    lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
    lammps_cells = {}
    for l in lammps_lines:
        l = filter(None, l)
        if (float(l[5]) + float(l[6]) + float(l[7])) > 0.0:
            lammps_cells[(float(l[1]), float(l[2]), float(l[3]))] = \
                        np.array([float(l[5]), float(l[6]), float(l[7])])

    if len(lammps_cells) != len(cfd_cells):
        print "Number of cells LAMMPS: ", len(lammps_cells), "\n",\
              "Number of cells dummy CFD: ", len(cfd_cells)
        print "FAILURE"
        sys.exit()

    for k in lammps_cells.keys():
        try:
            diff_forces = abs(cfd_cells[k] - lammps_cells[k])
            if (np.any(diff_forces > tol)):
                print "False"
                print cfd_cells[k]
                print lammps_cells[k]
                sys.exit()
        except KeyError:
            print "Cell not found."
            print k
            sys.exit()
    print "SUCCESS"

check_forces(1e-5)
