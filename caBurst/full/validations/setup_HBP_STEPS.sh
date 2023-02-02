#! /bin/bash

module load archive/2022-08

spack env create -d spackenv
sed -i '6 i\  concretization: together' spackenv/spack.yaml
spack env activate -d spackenv
spack add steps@develop^petsc/yjxaya
spack develop -p ${PWD} --no-clone steps@develop^petsc/yjxaya
spack install