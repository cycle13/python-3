#!/bin/sh
set -ex
F90=gfortran 

LIB_NETCDF="/usr/local/lib/"
INC_NETCDF="/usr/local/include/"


OMP='-lgomp'
F90FLAGS='-03'

cd ./src/fortran/

#PRE CLEAN
rm -f *.mod
rm -f *.o

#ln -sf ../../../common_python/common_functions/common_functions.f90        .
#ln -sf ../../../../../common/SFMT.f90          .

export FC=gfortran
export F90=gfortran

#f2py  -c --include-paths $INC_NETCDF -L$LIB_NETCDF -lnetcdff $OMP --f90flags="$F90FLAGS" -m scale_to_radar  common.f90 common_namelist_scale_to_radar.f90 scale_const.f90 scale_mapproj.F90 common_ncio.f90 common_scale.f90 common_obs_scale.f90 common_scale_to_radar.f90 

f2py  -c --include-paths $INC_NETCDF -L$LIB_NETCDF -lnetcdff $OMP --f90flags="$F90FLAGS" -m common_scale scale_const.f90 common_ncio.f90 common_scale.f90  


#CLEAN UP
rm -f *.mod
rm -f *.o
#rm common_functions.f90

echo "NORMAL END"
