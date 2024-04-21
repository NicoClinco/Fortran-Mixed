#!/bin/bash

g++ -c exampleCPP.cpp
gfortran -c exampleFOR.f90
gfortran -I. exampleCPP.o exampleFOR.o
