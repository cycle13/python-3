#!/bin/bash

id=$RANDOM

echo "#!/bin/bash            "        > run_slurm${id}.sh
echo "#SBATCH -p pps         "       >> run_slurm${id}.sh
echo "#SBATCH -J OpenMPI     "       >> run_slurm${id}.sh
echo "#SBATCH -t 24:00:00    "       >> run_slurm${id}.sh
echo "#SBATCH -c $2          "       >> run_slurm${id}.sh
echo "#SBATCH --mem 500000   "       >> run_slurm${id}.sh

#Very important! HDF5 read fails if this variable is set to true (the default!)
echo "HDF5_USE_FILE_LOCKING=\"FALSE\" " >> run_slurm${id}.sh
echo "mpiexec -np $2 python $1   "      >> run_slurm${id}.sh

sbatch run_slurm${id}.sh

