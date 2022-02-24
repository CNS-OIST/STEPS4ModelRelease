#! /bin/bash
# The line above is a "she-bang" to tell that /bin/bash is the command interpreting this file

#SBATCH --array=[1]
#SBATCH --partition=prod
#SBATCH --account=proj95
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2:00:00
#SBATCH --exclusive

set -x

module load archive/2021-12 python-dev
. ../../../oldspack/share/spack/setup-env.sh
spack load /saf25ux

nodes=$SLURM_JOB_NUM_NODES
ntasks=$(($nodes * 32))
seed=$(($SLURM_ARRAY_TASK_ID * 1))
time srun --nodes=$nodes --ntasks=$ntasks dplace python rallpack3.py $seed mesh/axon_cube_L1000um_D866nm_1135tets.msh 0

