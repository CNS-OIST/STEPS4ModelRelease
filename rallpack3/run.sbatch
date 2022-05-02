#! /bin/bash

# You need to set the appropriate SBATCH_PARTITION, SBATCH_ACCOUNT, and USE_STEPS4 mentioned in HBP_STEPS/doc/dev/README.md

#SBATCH --array=[1-10]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2:00:00
#SBATCH --exclusive
#SBATCH --mem=0

set -x

module load unstable python-dev python
spack load steps@develop

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))

time srun --nodes=$nodes --ntasks=$ntasks dplace \
python rallpack3.py $seed mesh/axon_cube_L1000um_D866nm_1135tets.msh $val_mon_steps_version

