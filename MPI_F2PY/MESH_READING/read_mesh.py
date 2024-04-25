import numpy as np
from mpi4py import MPI
#mpi4py.rc.initialize = False  # do not initialize MPI automatically
#mpi4py.rc.finalize = False    # do not finalize MPI automatically

#First direction of the research: The subroutines are
#called in the main program, while fortran is used to perform the loops, thus, it does not requires mpi, i think.
import global_vars

#fcomm = MPI.COMM_WORLD.py2f()
rank = MPI.COMM_WORLD.Get_rank()
#Read the mesh, suppose we have stored
#the data in a numpy array:
ncells = 10
nverts = 16
Cells = np.zeros((ncells,8),dtype=int)
Vertexes = np.zeros((nverts,3),dtype=int)
global_vars.global_variables.initialize_old_mesh(Vertexes,Cells,nverts,ncells)
global_vars.global_variables.initialize_new_mesh(Vertexes,Cells,nverts,ncells)

global_vars.global_variables.partitionmesh()

global_vars.global_variables.printnewlocalmesh()
