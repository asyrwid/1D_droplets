#!/bin/bash

# set parameters
Natoms="10"
L="52"
MaxOcc="10"
U="1"
r_vec="0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5"
#r_vec="0.05 0.075 0.1 0.125 0.15"
t="0.50"
MaxBondDim="512"
dir="./data/"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export ITENSOR_USE_OMP=1
#for t in $t_vec
#do
for r in $r_vec
	do
		nohup nice -19 ./droplet $Natoms $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_L."$L"_N."$Natoms"_U."$U"_t."$t"_r."$r &
#		./droplet $Natoms $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_t."$t"_r."$r &
done
#done
