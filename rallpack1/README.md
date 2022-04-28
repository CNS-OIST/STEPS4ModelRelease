# Rallpack 1: validations

In this model we can compare analytical, STEPS 3 and STEPS 4 solutions. `run.sbatch` is 
an example of a slurm job.

## Analytical solution

Since the model was extensively investigated in the past we hereby present the raw 
traces directly in: `results/analytical`                     
## STEPS 3

In order to run the STEPS 3 version of the model you just need to:

- Decide an rng seed. For example: 
`export rng_seed=1`
- Select a mesh file from the mesh folder. For example: 
`export mesh_file=mesh/axon_cube_L1000um_D866nm_1135tets.msh`
- Finally, run the command:
`python rallpack1_STEPS3.py ${rng_seed} ${mesh_file}`
- Raw traces will appear in the result folder: `results/STEPS3`
## STEPS 4

In order to run the STEPS 4 version of the model you just need to:

- Decide an rng seed. For example: 
`export rng_seed=1`
- Select a mesh file from the mesh folder. For example: 
`export mesh_file=mesh/axon_cube_L1000um_D866nm_1135tets.msh`
- Finally, run the command:
`python rallpack1_STEPS4.py ${rng_seed} ${mesh_file}`
- Raw traces will appear in the result folder: `results/STEPS4`