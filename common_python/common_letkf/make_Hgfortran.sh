#!/bin/bash

COMPILER=f2py
FFLAGS='-O3'
#FFLAGS='-O1 -fcheck=all'  #For debug

#This script compiles the fortran modules required to run the python experiments.

echo "Compiling DA routines"

rm -f *.mod *.o

$COMPILER -c -lgomp --opt="-fopenmp -lgomp" netlib.f90 SFMT.f90 common_tools.f90 common_mtx.f90 common_letkf.f90 common_da.f90 -m cletkf #> compile.out 2>&1

echo "Normal end"


