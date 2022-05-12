import numpy as np
import pickle
import sys
try:
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
except:
    "Error: PyFoam package is required to run the tests"

try:
    # Load parameters for the run
    params = pickle.load(open("md_params.dic", "rb"))

    # Parameters of the domain
    Lx = params["lx"]
    Ly = params["ly"]
    Lz = params["lz"]

    params = pickle.load(open("cfd_params.dic", "rb"))
    # Parameters of the mesh topology (cartesian grid)
    ncx = params["ncx"]
    ncy = params["ncy"]
    ncz = params["ncz"]
except:
    print("Error: Cannot read topology and domain data")


def compare_forces(tol, md_fname="md_forces.dat",
                   openfoam_dir="test_forces_case/2"):

    # OpenFOAM data   
    s_xy = ParsedParameterFile(openfoam_dir + "/sigmaxy")["internalField"]
    s_yz = ParsedParameterFile(openfoam_dir + "/sigmayz")["internalField"]
    s_xz = ParsedParameterFile(openfoam_dir + "/sigmaxz")["internalField"]
    s_xx = ParsedParameterFile(openfoam_dir + "/sigmaxx")["internalField"]
    s_yy = ParsedParameterFile(openfoam_dir + "/sigmayy")["internalField"]
    s_zz = ParsedParameterFile(openfoam_dir + "/sigmazz")["internalField"]
    cell_cx = ParsedParameterFile(openfoam_dir + "/ccx")["internalField"]
    cell_cy = ParsedParameterFile(openfoam_dir + "/ccy")["internalField"]
    cell_cz = ParsedParameterFile(openfoam_dir + "/ccz")["internalField"]

    openfoam_cells = {}
    cell_no = 0
    for cell_no in range(0, len(s_xy)): 
        cell_coord = (float(cell_cx[cell_no]), float(cell_cy[cell_no]), float(cell_cz[cell_no]))
        # Openfoam output cell centres with 6 decimal figures
        k = "{0:.5f}".format(cell_coord[0]), "{0:.5f}".format(cell_coord[1]),\
            "{0:.5f}".format(cell_coord[2])
        openfoam_cells[k] = np.array([float(s_xy[cell_no]), float(s_yy[cell_no]), 
                                      float(s_yz[cell_no])])


    # Line format of dummy md forces file -- > x y z sxy syy szy
    with open(md_fname, "r") as cfd_file:
        cfd_lines = cfd_file.readlines()
    md_lines = [l[:-1].split(" ") for l in cfd_lines]
    md_cells = {}
    for l in md_lines:
        k = "{0:.5f}".format(float(l[0])), "{0:.5f}".format(float(l[1])), "{0:.5f}".format(float(l[2]))
        md_cells[k] = np.array([float(l[3]), float(l[4]), float(l[5])])

    for k in list(md_cells.keys()):
        try:
            diff_forces = abs(md_cells[k] - openfoam_cells[k])
            if (np.any(diff_forces > tol)):
                print(md_cells[k])
                print(openfoam_cells[k])
                assert False
                sys.exit()
        except KeyError:
            print("Cell not found: cell " + str(k))
            assert False
            sys.exit()
    print("SUCCESS")

compare_forces(10e-6)



