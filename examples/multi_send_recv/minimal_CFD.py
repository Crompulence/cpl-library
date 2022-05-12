from mpi4py import MPI
from cplpy import CPL

comm = MPI.COMM_WORLD
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
cart_comm = CFD_COMM.Create_cart([1, 1, 1])
CPL.setup_cfd(cart_comm, xyzL=[1.0, 1.0, 1.0], 
              xyz_orig=[0.0, 0.0, 0.0], ncxyz=[32, 32, 32])
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=1)

for time in range(5):

    recv_array, ierr = CPL.recv(recv_array)
    print(("CFD", time, recv_array[0,0,0,0]))
    send_array[0,:,:,:] = 2.*time
    CPL.send(send_array)

CPL.finalize()

#Start again
CFD_COMM = CPL.init(CPL.CFD_REALM)
CPL.setup_cfd(cart_comm, xyzL=[1.0, 1.0, 1.0], 
              xyz_orig=[0.0, 0.0, 0.0], ncxyz=[32, 32, 32])
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=1)

for time in range(5):

    recv_array, ierr = CPL.recv(recv_array)
    print(("CFD", time, recv_array[0,0,0,0]))
    send_array[0,:,:,:] = 2.*time
    CPL.send(send_array)

CPL.finalize()

MPI.Finalize()
