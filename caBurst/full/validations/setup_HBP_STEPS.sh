#! /bin/bash

spack env deactivate

module purge

rm -rf spackenv

git stash
git checkout 0a7f75aa

spack env create -d spackenv
sed -i '6 i\  concretization: together' spackenv/spack.yaml
spack env activate -d spackenv
spack add steps@develop
spack develop -p ${PWD} --no-clone steps@develop
spack install