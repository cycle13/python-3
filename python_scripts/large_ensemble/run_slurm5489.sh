#!/bin/bash            
#SBATCH -p pps         
#SBATCH -J OpenMP      
#SBATCH -t 12:00:00    
#SBATCH -c 20          
#SBATCH --mem 30000    
source /opt/intel/bin/compilevars.sh intel64 
export OMP_NUM_THREADS=20 
python compute_update_comparison.py  
