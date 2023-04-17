#! /bin/bash

# You need to set the appropriate SBATCH_PARTITION, SBATCH_ACCOUNT, and USE_STEPS4 mentioned in HBP_STEPS/doc/dev/README.md

#SBATCH --array=[0-99%10]
#SBATCH --nodes=32
#SBATCH --time=00:30:00
#SBATCH --exclusive
#SBATCH --mem=0

set -x

module load unstable python-dev python
## load steps by hash
spack env activate -d /gpfs/bbp.cscs.ch/home/cattabia/projects/HBP_STEPS/spackenv
#export PYTHONPATH=`spack find --paths /$SPACK_INSTALLED_HASH | tail -n 1 | grep -o "/.*"`:$PYTHONPATH

steps_version=3

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))
mesh_fls=../../mesh/split_1024/steps3/CNG_segmented_2_split_1024.msh

time srun --nodes=$nodes --ntasks=$ntasks \
python caburst_short.py $seed $mesh_fls $steps_version


