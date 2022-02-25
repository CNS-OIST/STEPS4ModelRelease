#! /bin/bash

# You need to set the appropriate SBATCH_PARTITION, SBATCH_ACCOUNT mentioned in HBP_STEPS/doc/dev/README.md
# Sbatch does not accept $variables. $SBATCH_PARTITION and $SBATCH_ACCOUNT are here only for merely illustrative purposes.

#SBATCH --array=[1-100]
#SBATCH --partition=$SBATCH_PARTITION
#SBATCH --account=$SBATCH_ACCOUNT
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --exclusive

set -x

spack load steps

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))
mesh_fls=../../mesh/split_1024/steps3/CNG_segmented_2_split_1024.msh

time srun --nodes=$nodes --ntasks=$ntasks dplace python caBurstFullModel.py $seed $mesh_fls 0
