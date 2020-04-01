#!/bin/bash            
#SBATCH -J OpenMP      
#SBATCH -t 24:00:00    
#SBATCH -c 20          
#SBATCH --mem 30000    
export OMP_NUM_THREADS=20 
python compute_ensemble_statistics.py  
