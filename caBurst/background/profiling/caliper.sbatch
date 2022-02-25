#!/bin/bash

#SBATCH --account=
#SBATCH --partition=
#SBATCH --constraint=clx # Cascade Lake nodes on BB5 (intel codename clx)
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output="%x.out"

module purge
module load archive/2021-12
module load python-dev

# spack dev-build steps@develop+distmesh+petsc
export PYTHONPATH=~/STEPS/build/lib:$PYTHONPATH

spack load caliper # or module load caliper, depending on your setup
export CALI_CONFIG="runtime-report(calc.inclusive,output=stdout)"

export OMP_NUM_THREADS=1
for i in {1..30}
do
    echo " *** This is run $i *** "
    srun dplace python3 caBurstBackgroundModel.py 1 ../../mesh/split_${SLURM_NTASKS}/steps4/CNG_segmented_2_split_${SLURM_NTASKS} 1
    echo " *** This was run $i *** "
done
