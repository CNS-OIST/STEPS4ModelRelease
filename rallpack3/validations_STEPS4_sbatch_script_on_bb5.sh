#! /bin/bash

# You need to set the appropriate SBATCH_PARTITION, SBATCH_ACCOUNT mentioned in HBP_STEPS/doc/dev/README.md

#SBATCH --array=[1-1000]
#SBATCH --partition=$SBATCH_PARTITION
#SBATCH --account=$SBATCH_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2:00:00
#SBATCH --exclusive

set -x

spack load steps

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))

time srun --nodes=$nodes --ntasks=$ntasks dplace python rallpack3.py $seed mesh/axon_cube_L1000um_D866nm_1135tets.msh 1