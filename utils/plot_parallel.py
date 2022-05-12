import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import unittest


#def primefactors(num, f):

#    factors = np.zeros(num/2)

#    i = 2  #eligible factor
#    f = 1  #number of factors
#    n = num #store input number into a temporary variable
#    for
#        if (n%i == 0): #if i divides 2, it is a factor
#            factors[f] = i
#            f = f+1
#            n = n/i
#        else:
#            i = i+1     #not a factor. move to next number

#        if (n == 1):        
#            f = f-1        #its value will be one more than the number of factors
#            return


def collect_data(recv_array, cart_comm, ncxyz, plotrank=0):

    ncx, ncy, ncz = ncxyz
    rank = cart_comm.Get_rank()
    all_recv_array = cart_comm.gather(recv_array, root=plotrank)

    ncx_l = ncxyz[0]/cart_comm.Get_topo()[0][0]
    ncy_l = ncxyz[1]/cart_comm.Get_topo()[0][1]
    ncz_l = ncxyz[2]/cart_comm.Get_topo()[0][2]

    if rank == plotrank:
        field = np.zeros([recv_array.shape[0], ncx, ncy, ncz])
        #Loop over all processors
        for n, r in enumerate(all_recv_array):
            i, j, k = cart_comm.Get_coords(n)
            field[:, i*ncx_l:(i+1)*ncx_l, 
                     j*ncy_l:(j+1)*ncy_l, 
                     k*ncz_l:(k+1)*ncz_l] = r[:,:,:,:]

        return field
    else:
        return None


#Setup send and recv buffers
def allocate_buffer(ncxyz_l, cart_comm):
    ncx_l, ncy_l, ncz_l = ncxyz_l
    ir, jr, kr = cart_comm.Get_coords(cart_comm.Get_rank())
    recv_array = np.zeros((3, ncx_l, ncy_l, ncz_l), order='F', dtype=np.float64)
    for i in range(ncx_l):
        for j in range(ncy_l):
            for k in range(ncz_l):
                 recv_array[0,i,j,k] = i + ncx_l*ir
                 recv_array[1,i,j,k] = j + ncy_l*jr
                 recv_array[2,i,j,k] = k + ncz_l*kr

    return recv_array



if __name__ == "__main__":

    import unittest

    class TestPlot(unittest.TestCase):

        @classmethod
        def setUpClass(self):

            self.plotrank = 0

            #initialise MPI and CPL
            self.comm = MPI.COMM_WORLD

            # Parameters of the cpu topology (cartesian grid)
            ncx = 64; ncy = 64; ncz = 8
            self.ncxyz = [ncx, ncy, ncz]

            #primefactors(num, factors, f)

            self.npxyz = np.array([2, 2, 2], order='F', dtype=np.int32)
            self.xyzL = np.array([195.2503206, 18.62550553, 133.3416884], order='F', dtype=np.float64)
            self.xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

            #Setup coupled simulation
            self.cart_comm = self.comm.Create_cart([self.npxyz[0], self.npxyz[1], self.npxyz[2]])

            ncx_l = self.ncxyz[0]/self.npxyz[0]
            ncy_l = self.ncxyz[1]/self.npxyz[1]
            ncz_l = self.ncxyz[2]/self.npxyz[2]

            self.ncxyz_l = [ncx_l, ncy_l, ncz_l]

            self.recv_array = allocate_buffer(self.ncxyz_l, self.cart_comm)

        def test_collect(self):

            field = collect_data(self.recv_array, self.cart_comm, 
                                 self.ncxyz, self.plotrank)

            if self.cart_comm.Get_rank() == self.plotrank:
                for i in range(self.ncxyz[0]):
                    for j in range(self.ncxyz[1]):
                        for k in range(self.ncxyz[2]):
                            self.assertEqual(field[0,i,j,k], i) 
                            self.assertEqual(field[1,i,j,k], j)
                            self.assertEqual(field[2,i,j,k], k)

        def test_plot(self):

            field = collect_data(self.recv_array, self.cart_comm, 
                                 self.ncxyz, self.plotrank)

            if self.cart_comm.Get_rank() == self.plotrank:
                x = np.linspace(0.,1.,self.ncxyz[0])
                y = np.linspace(0.,1.,self.ncxyz[1])
                z = np.linspace(0.,1.,self.ncxyz[2])

                X, Y, Z = np.meshgrid(x, y, z)
                plt.pcolormesh(X[:,:,0], Y[:,:,0], 
                               np.mean(field[0,:,:,:],2), alpha=0.4)
                plt.colorbar()
                plt.ion()
                plt.show()
                plt.pause(2.)
                self.assertEqual(input("Plot looks correct? y/n:"),"y")
                plt.ioff()
        
        @classmethod
        def tearDownClass(self):
            self.cart_comm.Free()
            MPI.Finalize()
    

    unittest.main()

#    #initialise MPI and CPL
#    comm = MPI.COMM_WORLD
#    comm = comm
#    rank = comm.Get_rank()
#    nprocs_realm = comm.Get_size()

#    # Parameters of the cpu topology (cartesian grid)
#    ncx = 64; ncy = 64; ncz = 8
#    ncxyz = [ncx, ncy, ncz]
#    npxyz = np.array([2, 1, 2], order='F', dtype=np.int32)
#    NProcs = np.product(npxyz)
#    xyzL = np.array([195.2503206, 18.62550553, 133.3416884], order='F', dtype=np.float64)
#    xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)

#    #Setup coupled simulation
#    cart_comm = comm.Create_cart([npxyz[0], npxyz[1], npxyz[2]])

#    ncx_l = ncxyz[0]/npxyz[0]
#    ncy_l = ncxyz[1]/npxyz[1]
#    ncz_l = ncxyz[2]/npxyz[2]

#    recv_array = allocate_buffer(ncx_l, ncy_l, ncz_l, cart_comm)

#    plotrank = 0
#    field = collect_data(recv_array, cart_comm, ncxyz, plotrank)

#    if rank == plotrank:

#        for i in range(ncx):
#            for j in range(ncy):
#                for k in range(ncz):
#                    assert field[0,i,j,k] == i 
#                    assert field[1,i,j,k] == j
#                    assert field[2,i,j,k] == k


#        x = np.linspace(0.,1.,ncx)
#        y = np.linspace(0.,1.,ncy)
#        z = np.linspace(0.,1.,ncz)

#        X, Y, Z = np.meshgrid(x, y, z)
#        plt.pcolormesh(X[:,:,0], Y[:,:,0], np.mean(field[0,:,:,:],2), alpha=0.4)
#        plt.colorbar()
#        plt.show()

#        
#    MPI.Finalize()




