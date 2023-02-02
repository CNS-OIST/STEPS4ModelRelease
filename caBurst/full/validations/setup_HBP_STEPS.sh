#! /bin/bash

spack env deactivate

module purge

rm -rf spackenv



module load archive/2022-08

git checkout d846aba1bd212e03317cc376347a74582b82bc84

spack env create -d spackenv
sed -i '6 i\  concretization: together' spackenv/spack.yaml
spack env activate -d spackenv
spack add petsc@3.16.1
spack add steps@develop
spack develop -p ${PWD} --no-clone steps@develop^petsc@3.16.1
spack install