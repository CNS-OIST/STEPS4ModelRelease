# Rallpack 3: validations

In this model we can compare STEPS 3 and STEPS 4 solutions. It is suggested to run 
it on a cluster. `run.sbatch` is 
an example of a slurm job.
                 
## STEPS 3

In order to run the STEPS 3 version of the model you just need to:

- Decide an rng seed. For example: 
`export rng_seed=1`
- Select a mesh file from the mesh folder. For example: 
`export mesh_file=mesh/axon_cube_L1000um_D866nm_1135tets.msh`
- Finally, run the command:
`python rallpack1_STEPS3.py ${rng_seed} ${mesh_file}`
- Raw traces will appear in the result folder: `raw_traces/STEPS3`

## STEPS 4

In order to run the STEPS 4 version of the model you just need to:

- Decide an rng seed. For example: 
`export rng_seed=1`
- Select a mesh file from the mesh folder. For example: 
`export mesh_file=mesh/axon_cube_L1000um_D866nm_1135tets.msh`
- Finally, run the command:
`python rallpack1_STEPS4.py ${rng_seed} ${mesh_file}`
- Raw traces will appear in the result folder: `raw_traces/STEPS4`

