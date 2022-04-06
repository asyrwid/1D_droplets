#!/bin/bash

t_vec="$(seq 0.05 0.1 1)"
for t in $t_vec
do
 bash run_r_loop.sh $t
done
