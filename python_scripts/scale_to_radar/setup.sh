#!/bin/sh

#python3 setup.py build_ext --inplace
cd ./src/python/
f2py  -c -lgomp --f90flags="-fopenmp -lgomp -O3" -m calc_for ../fortran/calc_for.f90

f2py  -c -lgomp --f90flags="-fopenmp -lgomp -O3" -m radar_tools ../fortran/radar_tools.f90

