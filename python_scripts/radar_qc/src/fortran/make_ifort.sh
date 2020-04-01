#This script compiles the fortran code and generates a python module
export PATH=/opt/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=/opt/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=ifort
export F77=ifort
export F90=ifort
export F2PY=f2py3

#export FFLAGS='-fopenmp -lgomp -O3'
export FFLAGS='-g -traceback -fPIC'

$F2PY -c -lgomp --f90flags="$FFLAGS" -m common_functions common_functions.f90
$F2PY -c -lgomp --f90flags="$FFLAGS" -m common_qc_tools common_qc_tools.f90


