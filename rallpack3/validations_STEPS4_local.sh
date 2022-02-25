#! /bin/bash

set -x

nodes=1
tasks_per_node=32
seed=1
ntasks=$(($nodes * $tasks_per_node))

time srun --nodes=$nodes --ntasks=$ntasks dplace python rallpack3.py $seed mesh/axon_cube_L1000um_D866nm_1135tets.msh 1