#!/bin/bash

# set parameters
Na="10"
Nb="9"
L="60"
MaxOcc="10"
U="1"
r_vec="0.05 0.075 0.1 0.125 0.15 0.175"
#r_vec="0.2 0.225 0.25 0.275 0.3 0.325 0.35"
#r_vec="0.375 0.4 0.425 0.45 0.475 0.5"
t="1.00"
MaxBondDim="512"
epsilon="0.01"
dir="./data/"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export ITENSOR_USE_OMP=1

for r in $r_vec
	do
		nohup nice -19 ./droplet $Na $Nb $L $MaxOcc $t $U $r $MaxBondDim $dir $epsilon > "out_L."$L"_Na."$Na"_Nb."$Nb"_U."$U"_t."$t"_r."$r"_eps."$epsilon &
done

