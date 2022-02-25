#! /bin/bash

set -x

nodes=32
tasks_per_node=32
seed=1
ntasks=$(($nodes * $tasks_per_node))
mesh_fls=../../mesh/split_1024/steps4/CNG_segmented_2_split_1024

time srun --nodes=$nodes --ntasks=$ntasks dplace python caBurstFullModel.py $seed $mesh_fls 0
