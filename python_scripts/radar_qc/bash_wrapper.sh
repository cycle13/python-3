#!/bin/bash

my_script=$1

source activate gdal

ulimit -s unlimited

export OMP_NUM_THREADS=12 

python -u $my_script # 2>&1  ${my_script}.log
