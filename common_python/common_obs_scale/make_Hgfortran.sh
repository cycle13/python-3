#This script compiles the fortran code and generates a python module
export PATH=/share/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=/share/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=gfortran
export F90=gfortran

f2py  -c -lgomp --f90flags="-O3" -m common_obs_scale common_obs_scale.f90

