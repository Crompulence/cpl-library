from mpi4py import MPI
from cplpy import CPL

comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.MD_REALM)

nprocs = MD_COMM.Get_size()
rank = MD_COMM.Get_rank()

print(("MD code processor "+ str(rank+1) + " of " + str(nprocs)))

CPL.finalize()
MPI.Finalize()
