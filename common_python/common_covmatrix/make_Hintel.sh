#This script compiles the fortran code and generates a python module
export PATH=$HOME/data/software/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=$HOME/data/software/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=ifort
export F77=ifort
export F90=ifort

ln -sf ../common_functions/common_functions.f90      .
ln -sf ../common_random/common_random.f90            .
ln -sf ../common_random/SFMT.f90                     .

f2py  -c -lgomp --f90flags="-fopenmp -lgomp -O3" -m covariance_matrix_tools SFMT.f90 common_random.f90 common_functions.f90 covariance_matrix_tools.f90

#f2py3  -c --f90flags="-g" -m common_functions common_functions.f90

rm common_functions.f90
rm common_random.f90
rm SFMT.f90
