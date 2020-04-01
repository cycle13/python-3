#This script compiles the fortran code and generates a python module

#For paralel with openmp.
f2py -c -lgomp --f90flags="-fopenmp -lgomp" -m w_motion_vectors_old w_motion_vectors_parallel_old.f90

#Para ifort (esta fue  la compilacion del 21//9/2017)
#f2py/bin/f2py -c  -lgomp --f90flags="-fopenmp -lgomp -free" -m  w_motion_vectors w_motion_vectors_parallel.f90
#f2py/bin/f2py -c  -lgomp --f90flags="-fopenmp -lgomp" -m  w_motion_vectors w_motion_vectors_parallel.f90
#For single run without openmp.
#f2py -c -m motion_vectors motion_vectors_parallel.f90

