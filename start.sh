#!/bin/bash 

CutBondDimension_vec="512" #200 is minimal value
sweeps_tdvp_input_vec="4"

dmrg_N_max_sweeps="1000"
dmrg_nsweeps_vec="4"
on_off_time_evolution="0"
#parameters for (MPO,MPO) multiplication
nmultMPO_cutoff_vec="10" #10^{-#}
nmultMPO_MaxDim_vec="500" #before was 500
nmultMPO_method_vec="DensityMatrix"
epsilon_K_tdvp_vec="10"
add_basis_method_tdvp_vec="DensityMatrix"
NoOfSteps_tdvp_vec="5"


PBC_vec="1" # 1: Periodic Boundary Conditions; 0: Open Boundary Conditions


U_vec="1"
L_vec="20"
N_particles_vec="6 8 4 10 16"
N_particles_vec="4"


t_vec="1"
r_vec="0.10"


#Imaginary tdvp time evolution parameters
N_max_tdvp="10"
dt_tdvp="0.01"

#Real time evolution parameters
t_max_vec="100"
N_time_shot_vec="100"


 

 
for NoOfSteps_tdvp in $NoOfSteps_tdvp_vec
do
for PBC in $PBC_vec
do
for L in $L_vec
do                                           
	echo L=$L
for N_particles in $N_particles_vec
do
	echo N=$N_particles
for t in $t_vec
do	
for U in $U_vec
do
for r in $r_vec
do                               	
for dmrg_nsweeps in $dmrg_nsweeps_vec
do                 
for CutBondDimension in $CutBondDimension_vec
do
for sweeps_tdvp_input in $sweeps_tdvp_input_vec
do
for nmultMPO_cutoff in $nmultMPO_cutoff_vec
do
for nmultMPO_MaxDim in $nmultMPO_MaxDim_vec
do
for nmultMPO_method in $nmultMPO_method_vec
do
for epsilon_K_tdvp in $epsilon_K_tdvp_vec
do
for add_basis_method_tdvp in $add_basis_method_tdvp_vec
do
((maxOccupation=$N_particles/2))
echo $maxOccupation
	time ./bose_bose_mixture $CutBondDimension $sweeps_tdvp_input  $nmultMPO_cutoff $nmultMPO_MaxDim $nmultMPO_method $epsilon_K_tdvp $add_basis_method_tdvp $PBC $L $N_particles $t $U $r  $maxOccupation $NoOfSteps_tdvp $dmrg_N_max_sweeps
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done






