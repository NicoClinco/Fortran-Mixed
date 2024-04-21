#!/bin/bash

f2py --f90exec=mpifort --f77exec=mpifort -c subroutine.f90 -m simpleMPI


