#! /bin/bash
# The line above is a "she-bang" to tell that /bin/bash is the command interpreting this file

#SBATCH --array=[1-100]
#SBATCH --partition=prod
#SBATCH --account=proj95
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --exclusive

set -x

module load archive/2021-12 python-dev
. /gpfs/bbp.cscs.ch/project/proj16/cattabia/oldspack/share/spack/setup-env.sh
spack load /saf25ux


for nodes in ${SLURM_JOB_NUM_NODES}; do
    ntasks=$(($nodes * 32))
    seed=$(($SLURM_ARRAY_TASK_ID * 1))
    MESH_FLS=/gpfs/bbp.cscs.ch/project/proj16/cattabia/projects/STEPS4Models/CaBurstUnifiedMesh/split_1024/steps3/CNG_segmented_2_split_1024.msh
    time srun --nodes=$nodes --ntasks=$ntasks dplace \
 python ../caburstfull_single_script.py ${seed} $MESH_FLS 0
done

