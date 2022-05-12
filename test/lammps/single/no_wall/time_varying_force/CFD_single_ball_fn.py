import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpi4py import MPI

from cplpy import CPL


def CFD(xyzL = [1.5E-003, 1.5E-003, 2.50E-003], 
        g = 9.81, ncxyz=[8,8,8], npxyz=[1,1,1],
        Nsteps=101):

    #initialise MPI and CPL
    comm = MPI.COMM_WORLD
    CPL = CPL()
    MD_COMM = CPL.init(CPL.CFD_REALM)
    nprocs_realm = MD_COMM.Get_size()

    ## Parameters of the cpu topology (cartesian grid)
    npxyz = np.array(npxyz, order='F', dtype=np.int32)
    NProcs = np.product(npxyz)

    xyzL = np.array(xyz, order='F', dtype=np.float64)
    xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
    ncxyz = np.array(ncxyz, order='F', dtype=np.int32)
    if (nprocs_realm != NProcs):
        print(("Non-coherent number of processes in MD ", nprocs_realm,
                " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2]))
        MPI.Abort(errorcode=1)

    #Setup coupled simulation
    cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
    CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

    #Get constraint region
    cnst_limits = CPL.get_cnst_limits();
    cnst_portion = CPL.my_proc_portion(cnst_limits)
    [cnst_ncxl, cnst_ncyl, cnst_nczl] = CPL.get_no_cells(cnst_portion)

    #Get overlap region
    olap_limits = CPL.get_olap_limits()
    BC_limits = np.array([olap_limits[0], olap_limits[1], 
                          olap_limits[2], olap_limits[3], 
                          olap_limits[4], olap_limits[5]], dtype=np.int32)
    BC_portion = CPL.my_proc_portion(BC_limits)
    [BC_ncxl, BC_ncyl, BC_nczl] = CPL.get_no_cells(BC_portion)

    #Allocate send and recv arrays
    recv_array = np.zeros((4, BC_ncxl, BC_ncyl, BC_nczl), order='F', dtype=np.float64)
    send_array = np.zeros((9, cnst_ncxl, cnst_ncyl, cnst_nczl), order='F', dtype=np.float64)

    for time in range(Nsteps):

        # send data to update
        send_array[2,:,:,:] = -5.9490638385009208e-08*g  # * mi
        CPL.send(send_array, cnst_portion)

        # recv data and plot
        recv_array, ierr = CPL.recv(recv_array, BC_portion)

        print(time)

    CPL.finalize()
    MPI.Finalize()




