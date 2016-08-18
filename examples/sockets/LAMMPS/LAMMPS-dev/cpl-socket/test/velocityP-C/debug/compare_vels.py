import sys
import numpy as np


def compare_vels(tol, lammps_fname="lammps_vels.dat", cfd_fname="cfd_vels.dat"):

    # Line format of CFD script file -- > x z vx vy vz
    with open(cfd_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
    cfd_cells = {}
    for l in cfd_lines:
        cfd_cells[(float(l[0]), float(l[1]))] = np.array([float(l[2]),
                                                         float(l[3]),
                                                         float(l[4])])

    # Line format of LAMMPS file -- > chunk x z ncount vx vy vz
    with open(lammps_fname, "r") as lammps_file:
        lammps_lines = lammps_file.readlines()
    skip = int(lammps_lines[3].split(" ")[1])
    lammps_lines = lammps_lines[4:]
    lammps_lines = lammps_lines[:skip]
    lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
    lammps_cells = {}
    for l in lammps_lines:
        l = filter(None, l)
        lammps_cells[(float(l[1]), float(l[2]))] = np.array([float(l[4]),
                                                            float(l[5]),
                                                            float(l[6])])
    if len(lammps_cells) != len(cfd_cells):
        print "FALSE"
        sys.exit()

    for k in lammps_cells.keys():
        try:
            diff_vel = abs(cfd_cells[k] - lammps_cells[k])
            if (np.any(diff_vel > tol)):
                print "False"
                print cfd_cells[k]
                print lammps_cells[k]
                sys.exit()
        except KeyError:
            print "Cell not found."
            print k
            sys.exit()
    print "SUCCESS"
compare_vels(1e-6)
