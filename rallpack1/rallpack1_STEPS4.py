# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack3 model
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import math
import os
import pandas as pd
import sys

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# MESHFILE = 'axon_cube_L1000um_D866nm_1135tets.msh'
MESHFILE = sys.argv[2]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Leak conductance, Siemens/m^2
L_G = 0.25

# Leak reveral potential, V
leak_rev = -65.0e-3

# Total leak conductance for ideal cylinder:
surfarea_cyl = 1.0*math.pi*1000*1e-12
L_G_tot = L_G*surfarea_cyl

# Ohm.m
Ra = 1.0

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # 

# Sim end time (seconds)
SIM_END = 0.25

# The current injection in amps
Iinj = 0.1e-9

SAVE_DT = 5e-5
# EFIELD DT
EF_DT = 1e-6

SEED = int(sys.argv[1])

# # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # # 

# record potential at the two extremes along (z) axis 
POT_POS = [0.0, 1.0e-03]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Mesh geometry
mesh = DistMesh(MESHFILE)

with mesh:
    cyto = Compartment.Create()

    memb = Patch.Create(cyto, None, 'ssys')
    z_min = Patch.Create(cyto, None, 'ssys')
    z_max = Patch.Create(cyto, None, 'ssys')

    surfarea_mesh = memb.Area
    corr_fac_area = surfarea_mesh/surfarea_cyl

    vol_cyl = math.pi*0.5*0.5*1000*1e-18
    vol_mesh = cyto.Vol
    corr_fac_vol = vol_mesh/vol_cyl

    membrane = Membrane.Create([memb], capacitance=0.01 / corr_fac_area)

    cyto.Conductivity = 1 / (Ra * corr_fac_vol)
    
    # The tetrahedrons from which to record potential
    POT_TET = TetList(mesh.tets[0, 0, z] for z in POT_POS)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mdl = Model()
r = ReactionManager()
with mdl:
    ssys = SurfaceSystem.Create()

    # Leak
    leaksus = SubUnitState.Create()
    Leak = Channel.Create([leaksus])

    with ssys:
        # Set the single-channel conductance:
        g_leak_sc = L_G_tot/len(membrane.tris)
        OC_L = OhmicCurr.Create(Leak[leaksus], g_leak_sc, leak_rev) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Create the solver objects
rng = RNG('mt19937', 512, SEED)
sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)

# Data saving
rs = ResultSelector(sim)

Vrs = rs.TETS(POT_TET).V

sim.toSave(Vrs, dt=SAVE_DT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim.newRun()

sim.TRIS(membrane.tris).Leak[leaksus].Count = 1

sim.membrane.Potential = -65e-3

minzverts = list(set([v for t in z_min.tris for v in t.verts]))
for v in minzverts:
    sim.solver.setVertIClamp(v.idx, Iinj/len(minzverts))

# sim.stepsSolver.setEfieldTolerances(1e-50, 1e-8)

sim.EfieldDT = EF_DT

progress_dt = SAVE_DT*10
nsteps = math.ceil(SIM_END/progress_dt)
for i in range(nsteps):
    current_time = i*progress_dt
    sim.run(current_time)
    print(f"Progress: {round(1e4*(i/nsteps))/1e2}%")

if MPI.rank == 0:
    df = pd.DataFrame({"t":Vrs.time[0], "V_z_min":Vrs.data[0,:,0], "_Vz_max":Vrs.data[0,:,1]})
    df.to_csv(f'results/STEPS4/res{SEED}_STEPS4_{len(mesh.verts)*3}DoFs.txt', sep=" ", index=False)

