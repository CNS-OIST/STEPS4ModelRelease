#! /bin/bash

# You need to set the appropriate SBATCH_PARTITION, SBATCH_ACCOUNT, and USE_STEPS4 mentioned in HBP_STEPS/doc/dev/README.md

#SBATCH --array=[1-100]
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=64
#SBATCH --time=2:00:00
#SBATCH --exclusive
#SBATCH --mem=0

set -x

module load unstable python-dev python
#spack load steps@develop
## load steps by hash
export PYTHONPATH=`spack find --paths /$SPACK_INSTALLED_HASH | tail -n 1 | grep -o "/.*"`:$PYTHONPATH

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))
mesh_fls=../../mesh/split_1024/steps3/CNG_segmented_2_split_1024.msh

time srun --nodes=$nodes --ntasks=$ntasks dplace \
python caBurstFullModel.py $seed $mesh_fls $steps_version

