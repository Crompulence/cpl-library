import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpi4py import MPI

from cplpy import CPL

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)
nprocs_realm = MD_COMM.Get_size()



## Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
NProcs = np.product(npxyz)

xyzL = np.array([1.5000000000000000E-003, 1.5000000000000000E-003, 2.5000000000000001E-003], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([8, 8, 8], order='F', dtype=np.int32)
if (nprocs_realm != NProcs):
    print("Non-coherent number of processes in MD ", nprocs_realm,
            " no equal to ",  npxyz[0], " X ", npxyz[1], " X ", npxyz[2])
    MPI.Abort(errorcode=1)

#Setup coupled simulation
cart_comm = MD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

##Plot output
fig, axs = plt.subplots(2,1)
axs[0].set_ylabel('$u$')
axs[1].set_ylabel(r'$\rho$')
#plt.subplots_adjust(right=0.9)
#axslider = plt.axes([0.25, 0.1, 0.65, 0.03])
freq = 0
#freq = 0.25; vmin = -2.0; vmax = 2.0
#sfreq = Slider(axslider, '$U_{wall}$', vmin, vmax, valinit=freq)
#def update(val):
#    freq = sfreq.val
#    global freq
#    print("CHANGED", freq)
#sfreq.on_changed(update)

plt.ion()
plt.show()

## === Plot both grids ===
dx = CPL.get("xl_md")/float(CPL.get("ncx"))
dy = CPL.get("yl_md")/float(CPL.get("ncy"))
dz = CPL.get("zl_md")/float(CPL.get("ncz"))
dV = dx*dy*dz

#first_time = True


cnst_limits = CPL.get_cnst_limits();
cnst_portion = CPL.my_proc_portion(cnst_limits)
[cnst_ncxl, cnst_ncyl, cnst_nczl] = CPL.get_no_cells(cnst_portion)

olap_limits = CPL.get_olap_limits()
BC_limits = np.array([olap_limits[0], olap_limits[1], 
                      olap_limits[2], olap_limits[3], 
                      olap_limits[4], olap_limits[5]], dtype=np.int32)

BC_portion = CPL.my_proc_portion(BC_limits)
[BC_ncxl, BC_ncyl, BC_nczl] = CPL.get_no_cells(BC_portion)
recv_array = np.zeros((4, BC_ncxl, BC_ncyl, BC_nczl), order='F', dtype=np.float64)

ft = True
for time in range(100000):

    # send data to update
    send_array = freq*np.ones((9, cnst_ncxl, cnst_ncyl, cnst_nczl), order='F', dtype=np.float64)
    send_array[2,:,:,:] = -5.9490638385009208e-08*9.81# * mi
    CPL.send(send_array, cnst_portion)

    # recv data and plot
    recv_array, ierr = CPL.recv(recv_array, BC_portion)
    print("Python recvs", time, recv_array.shape, np.min(recv_array), np.max(recv_array), np.sum(recv_array))

    ##allcplpass##################
#    for i in range(recv_array.shape[1]):
#        for k in range(recv_array.shape[3]):
#            print(time,i,k,recv_array[0,i,:,k],recv_array[3,i,:,k],np.mean(recv_array[0,:,:,:],(0,2))/(dV*0.6))
#            l = ax.plot(recv_array[0,i,:,k]/(dV*0.5),'-o')
#    l = ax.plot(recv_array[2,:,0,:],recv_array[2,:,0,:],'-o',lw=3.)
#    ax.set_ylim([-2.,2.])
    l = axs[0].pcolormesh(np.mean(recv_array[2,:,:,:],1).T)
    l = axs[1].pcolormesh(np.mean(recv_array[3,:,:,:],1).T)

#    if ft:
#        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#        fig.colorbar(l, cax=cax)
#        ft=False
    ##allcplpass##################
    plt.pause(0.001)
    [ax.cla() for ax in axs]



CPL.finalize()
MPI.Finalize()




