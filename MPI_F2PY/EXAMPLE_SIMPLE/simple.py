import mpi4py
mpi4py.rc.initialize = False  # do not initialize MPI automatically
mpi4py.rc.finalize = False    # do not finalize MPI automatically

#from mpi4py import MPI
import simpleMPI

simpleMPI.simple()



