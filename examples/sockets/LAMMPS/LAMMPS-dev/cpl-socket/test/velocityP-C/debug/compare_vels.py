import sys
import numpy as np


def compare_vels(tol, lammps_fname="lammps_vels.dat", cfd_fname="cfd_vels.dat"):

    # Line format of CFD script file -- > x y z vx vy vz
    with open(cfd_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
    cfd_cells = {}
    for l in cfd_lines:
        cfd_cells[(float(l[0]), float(l[1]), float(l[2]))] = np.array([float(l[3]),
                                                                       float(l[4]),
                                                                       float(l[5])])

    # Line format of LAMMPS file -- > chunk x y z ncount vx vy vz
    with open(lammps_fname, "r") as lammps_file:
        lammps_lines = lammps_file.readlines()
    skip = int(lammps_lines[3].split(" ")[1])
    lammps_lines = lammps_lines[4:]
    lammps_lines = lammps_lines[:skip]
    lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
    lammps_cells = {}
    for l in lammps_lines:
        l = filter(None, l)
        lammps_cells[(float(l[1]), float(l[2]), float(l[3]))] = np.array([float(l[5]),
                                                                          float(l[6]),
                                                                          float(l[7])])
    for k in cfd_cells.keys():
        try:
            diff_vel = abs(cfd_cells[k] - lammps_cells[k])
            if (np.any(diff_vel > tol)):
                print "Cell mismatch for key: " + str(k)
                print "CFD cell: " + str(cfd_cells[k])
                print "LAMMPS cell: " + str(lammps_cells[k])
                sys.exit()
        except KeyError:
            print "Cell not found :" + str(k)
            sys.exit()
    print "SUCCESS"
compare_vels(1e-6)
