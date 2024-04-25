import mpi4py
mpi4py.rc.initialize = False  # do not initialize MPI automatically
mpi4py.rc.finalize = False    # do not finalize MPI automatically

from mpi4py import MPI 

MPI.Init()

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

print(f"We are on rank: {my_rank}")
MPI.Finalize()  


