#!/bin/bash

# set parameters
Natoms="10"
L="102"
MaxOcc="5"
U="1"
r_vec="$(seq 0.05 0.025 0.5)"
MaxBondDim="256"
dir="./data/"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export ITENSOR_USE_OMP=1
t=$1 #take a passed argument from run_t_loop.sh script
for r in $r_vec
do
 nohup nice -19 ./droplet $Natoms $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_t."$t"_r."$r &
done
