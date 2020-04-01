#!/bin/bash            
#SBATCH -p pps         
#SBATCH -J OpenMPI     
#SBATCH -t 24:00:00    
#SBATCH -c 50          
#SBATCH --mem 500000   
HDF5_USE_FILE_LOCKING="FALSE" 
mpiexec -np 50 python convert_letkfout_LE_D1_1km_30sec.py   
