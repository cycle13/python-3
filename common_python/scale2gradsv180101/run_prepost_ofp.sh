#!/bin/bash            
source activate basemap
#mpiexec -np 68 python ./convert_letkfout_LE_D1_1km_30sec.py  
mpiexec -np 1 python ./convert_letkfout_LE_D1_1km_2min.py 
