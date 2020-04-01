#!/bin/bash            
#SBATCH -p pps         
#SBATCH -J OpenMPI     
#SBATCH -t 24:00:00    
#SBATCH -c 30          
#SBATCH --mem 500000   
mpiexec -np 30 python convert_letkfout_LE_D1_1km_30sec_nospinup.py   
