#!/bin/bash

id=$RANDOM

echo "#!/bin/bash            "        > run_slurm${id}.sh
echo "#SBATCH -p pps         "       >> run_slurm${id}.sh
echo "#SBATCH -J OpenMP      "       >> run_slurm${id}.sh
echo "#SBATCH -t 24:00:00    "       >> run_slurm${id}.sh
echo "#SBATCH -c 20          "       >> run_slurm${id}.sh
echo "#SBATCH --mem 30000    "       >> run_slurm${id}.sh

#echo "source /opt/intel/bin/compilevars.sh intel64 "  >> run_slurm${id}.sh

echo "export OMP_NUM_THREADS=$2 "    >> run_slurm${id}.sh

echo "python $1  "                   >> run_slurm${id}.sh


sbatch run_slurm${id}.sh      

