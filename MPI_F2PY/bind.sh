#!/bin/bash

f2py --f90exec=mpif90 --f77exec=mpif90 -c subroutine.f90 -m simpleMPI \
     -L/usr/lib/x86_64-linux-gnu/ -lmpi
#f2py --f90exec=mpif90 --f77exec=mpif90 -I
#mpiexec -n 2 python3 mpi4py_example.py

