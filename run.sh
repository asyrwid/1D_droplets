#!/bin/bash

# set parameters
Natoms="1"
L="10"
MaxOcc="2"
U="1"
#r_vec="0.1 0.15 0.2 0.25"
#t_vec="0.05 0.1 0.2 0.3"
r_vec="0.2"
t_vec="0.5"
MaxBondDim="256"
dir="./data/"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export ITENSOR_USE_OMP=1
for t in $t_vec
do
for r in $r_vec
	do
	#	nohup nice -19 ./droplet $Natoms $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_L."$L"_N."$Natoms"_U."$U"_t."$t"_r."$r &

		./droplet $Natoms $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_t."$t"_r."$r &
done
done
