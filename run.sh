#!/bin/bash

# set parameters
Natoms="2"
L="50"
MaxOcc="2"
t="1."
U="5."
r="0.1 0.15 0.2 0.25"
MaxBondDim="256"
dir="/home/asyrwid/ITensor-3.1.6/programs/1D_droplets/data/"

for R in $r
	do
		nohup nice -19 ./droplet $Natoms $L $MaxOcc $t $U $R $MaxBondDim $dir > "out_r="$R &
done
