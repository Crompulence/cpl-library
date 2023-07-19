import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

from cplpy import CPL
from draw_grid import draw_grid


class CFD():


    def __init__(self, npxyz, xyzL, xyz_orig, ncxyz):

        #initialise MPI and CPL
        self.comm = MPI.COMM_WORLD
        self.CPL = CPL()
        self.CFD_COMM = self.CPL.init(CPL.CFD_REALM)
        self.nprocs_realm = self.CFD_COMM.Get_size()

        # Parameters of the cpu topology (cartesian grid)
        self.npxyz = np.array(npxyz, order='F', dtype=np.int32)
        self.NProcs = np.product(npxyz)
        self.xyzL = np.array(xyzL, order='F', dtype=np.float64)
        self.xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)
        self.ncxyz = np.array(ncxyz, order='F', dtype=np.int32)

        if (self.nprocs_realm != self.NProcs):
            print(("Non-coherent number of processes in CFD ", self.nprocs_realm,
                    " no equal to ",  self.npxyz[0], " X ", self.npxyz[1], " X ", self.npxyz[2]))
            MPI.Abort(errorcode=1)

        #Setup coupled simulation
        self.cart_comm = self.CFD_COMM.Create_cart([self.npxyz[0], self.npxyz[1], self.npxyz[2]])
        self.CPL.setup_cfd(self.cart_comm, self.xyzL, self.xyz_orig, self.ncxyz)

        #Get limits of overlap region
        self.olap_limits = self.CPL.get_olap_limits()
        self.portion = self.CPL.my_proc_portion(self.olap_limits)
        [self.ncxl, self.ncyl, self.nczl] = self.CPL.get_no_cells(self.portion)

        self.dx = self.CPL.get("xl_cfd")/float(self.CPL.get("ncx"))
        self.dy = self.CPL.get("yl_cfd")/float(self.CPL.get("ncy"))
        self.dz = self.CPL.get("zl_cfd")/float(self.CPL.get("ncz"))
        self.ioverlap = (self.CPL.get("icmax_olap")-self.CPL.get("icmin_olap")+1)
        self.joverlap = (self.CPL.get("jcmax_olap")-self.CPL.get("jcmin_olap")+1)
        self.koverlap = (self.CPL.get("kcmax_olap")-self.CPL.get("kcmin_olap")+1)
        self.xoverlap = self.ioverlap*self.dx
        self.yoverlap = self.joverlap*self.dy
        self.zoverlap = self.koverlap*self.dz

    def recv_CPL_data(self):

        # recv data to plot
        self.recv_array = np.zeros((1, self.ncxl, self.ncyl, self.nczl), order='F', dtype=np.float64)
        self.recv_array, ierr = self.CPL.recv(self.recv_array, self.olap_limits)

    def plot_grid(self, ax):

        #Plot CFD and coupler Grid
        draw_grid(ax, 
                  nx=self.CPL.get("ncx"),
                  ny=self.CPL.get("ncy"),
                  nz=self.CPL.get("ncz"),
                  px=self.CPL.get("npx_cfd"),
                  py=self.CPL.get("npy_cfd"),
                  pz=self.CPL.get("npz_cfd"),
                  xmin=self.CPL.get("x_orig_cfd"),
                  ymin=self.CPL.get("y_orig_cfd"),
                  zmin=self.CPL.get("z_orig_cfd"),
                  xmax=(self.CPL.get("icmax_olap")+1)*self.dx,
                  ymax=self.CPL.get("yl_cfd"),
                  zmax=(self.CPL.get("kcmax_olap")+1)*self.dz,
                  lc = 'r',
                  label='CFD')

        #Plot MD domain
        draw_grid(ax, nx=1, ny=1, nz=1,
                  px=self.CPL.get("npx_md"),
                  py=self.CPL.get("npy_md"),
                  pz=self.CPL.get("npz_md"),
                  xmin=self.CPL.get("x_orig_md"),
                  ymin=-self.CPL.get("yl_md")+self.yoverlap,
                  zmin=self.CPL.get("z_orig_md"),
                  xmax=(self.CPL.get("icmax_olap")+1)*self.dx,
                  ymax=self.yoverlap,
                  zmax=(self.CPL.get("kcmax_olap")+1)*self.dz,
                  label='MD')


    def plot_data(self, ax):

        # === Plot both grids ===

        #Plot x component on grid
        x = np.linspace( self.portion[0]*self.dx+.5*self.dx,
                        (self.portion[0]+self.ncxl)*self.dx-.5*self.dx,self.ncxl)
        z = np.linspace(self.portion[0]*self.dx+.5*self.dz,
                        self.nczl*self.dx-.5*self.dz,self.nczl)


        try:
            for j in range(self.joverlap):
                ax.plot(x, 0.5*self.dy*(self.recv_array[0,:,j,0]+1.+2*j), 's-')
        except ValueError:
            print(("Arrays not equal:", x.shape, z.shape, self.recv_array.shape))
            raise

    def finalise(self):

        self.CPL.finalize()
        MPI.Finalize()


if __name__ == '__main__':

    #Get input file
    import inpututils
    ip = inpututils.InputMod("./CFD.in")
    npxyz = ip.read_input("npxyz")
    xyzL = ip.read_input("xyzL")
    xyz_orig = ip.read_input("xyz_orig")
    ncxyz = ip.read_input("ncxyz")

    cfd = CFD(npxyz=npxyz, 
              xyzL = xyzL, 
              xyz_orig = xyz_orig, 
              ncxyz = ncxyz)
    cfd.recv_CPL_data()
    fig, ax = plt.subplots(1,1)
    cfd.plot_grid(ax)
    cfd.plot_data(ax)
    plt.show()
    cfd.finalise()
    

