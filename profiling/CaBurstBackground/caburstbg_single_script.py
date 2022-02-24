########################################################################                                                                          
import steps.interface
                                                                                                                                                  
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

from extra.constants import *

import numpy as np
import pandas as pd
import sys

import psutil
mem_ini = psutil.Process().memory_info().rss / 1024**2

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if MPI.rank == 0:
    print("*-*-*->", "Import : Done", flush=True)

SEED = int(sys.argv[1])
mesh_file = sys.argv[2]
USE_STEPS_4 = int(sys.argv[3]) > 0

SIM_TIME = 30.0e-5

if MPI.rank == 0:
    if USE_STEPS_4:
        print("*-*-*->", "STEPS4 SIMULATION", flush=True)
    else:
        print("*-*-*->", "STEPS3 SIMULATION", flush=True)

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
r = ReactionManager()
with mdl:
    # Species
    Ca, Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

    # Vol/surface systems
    vsys = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()

    with vsys:
        diff_Ca =     Diffusion.Create(Ca, DCST)
        diff_CBsf =   Diffusion.Create(CBsf, DCB)
        diff_CBsCa =  Diffusion.Create(CBsCa, DCB)
        diff_CBCaf =  Diffusion.Create(CBCaf, DCB)
        diff_CBCaCa = Diffusion.Create(CBCaCa, DCB)
        diff_PV =     Diffusion.Create(PV, DPV)
        diff_PVCa =   Diffusion.Create(PVCa, DPV)
        diff_PVMg =   Diffusion.Create(PVMg, DPV)

        # iCBsf fast and slow
        (iCBsf + Ca <r[1]> iCBsCa) + Ca <r[2]> iCBCaCa
        (iCBsf + Ca <r[3]> iCBCaf) + Ca <r[4]> iCBCaCa
        r[1].K = iCBsf1_f_kcst, iCBsf1_b_kcst
        r[2].K = iCBsCa_f_kcst, iCBsCa_b_kcst
        r[3].K = iCBsf2_f_kcst, iCBsf2_b_kcst
        r[4].K = iCBCaf_f_kcst, iCBCaf_b_kcst

        # CBsf fast and slow
        (CBsf + Ca <r[1]> CBsCa) + Ca <r[2]> CBCaCa
        (CBsf + Ca <r[3]> CBCaf) + Ca <r[4]> CBCaCa
        r[1].K = CBsf1_f_kcst, CBsf1_b_kcst
        r[2].K = CBsCa_f_kcst, CBsCa_b_kcst
        r[3].K = CBsf2_f_kcst, CBsf2_b_kcst
        r[4].K = CBCaf_f_kcst, CBCaf_b_kcst

        # PVCa
        PV + Ca <r[1]> PVCa
        r[1].K = PVca_f_kcst, PVca_b_kcst

        # PVMg
        PV + Mg <r[1]> PVMg
        r[1].K = PVmg_f_kcst, PVmg_b_kcst

        # Ca Influx converted from P Type current
        None >r['CaInflux']> Ca
        r['CaInflux'].K = 0.0

    with ssys:
        # Ca Pump
        Pump.s + Ca.i <r[1]> CaPump.s >r[2]> Pump.s
        r[1].K = P_f_kcst, P_b_kcst
        r[2].K = P_k_kcst
    
########### MESH & COMPARTMENTALIZATION #################

mesh = DistMesh(mesh_file, 1e-6) if USE_STEPS_4 else TetMesh.LoadGmsh(mesh_file, 1e-6)

with mesh:
    if USE_STEPS_4:
        __MESH__ = Compartment.Create(vsys)

        __MESH_BOUNDARY__ = Patch.Create(__MESH__, None, ssys)
    else:
        __MESH__ = Compartment.Create(mesh.tets, vsys)

        __MESH_BOUNDARY__ = Patch.Create(mesh.surface, __MESH__, None, ssys)

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, SEED)

if USE_STEPS_4:
    sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
else:
    part = GmshPartition(mesh)
    sim = Simulation('TetOpSplit', mdl, mesh, rng, MPI.EF_NONE, part)

sim.newRun()

surfarea = __MESH_BOUNDARY__.Area
pumpnbs = 6.022141e12*surfarea

sim.__MESH_BOUNDARY__.Pump.Count = round(pumpnbs)
sim.__MESH_BOUNDARY__.CaPump.Count = 0

sim.__MESH__.Mg.Conc = Mg_conc

sim.__MESH__.iCBsf.Conc = iCBsf_conc
sim.__MESH__.iCBCaf.Conc = iCBCaf_conc
sim.__MESH__.iCBsCa.Conc = iCBsCa_conc
sim.__MESH__.iCBCaCa.Conc = iCBCaCa_conc

sim.__MESH__.CBsf.Conc = CBsf_conc
sim.__MESH__.CBCaf.Conc = CBCaf_conc
sim.__MESH__.CBsCa.Conc = CBsCa_conc
sim.__MESH__.CBCaCa.Conc = CBCaCa_conc

sim.__MESH__.PV.Conc = PV_conc
sim.__MESH__.PVCa.Conc = PVCa_conc
sim.__MESH__.PVMg.Conc = PVMg_conc

if MPI.rank == 0:
    print("Simulating model, it will take a while if running with small amount of processes...")

import time
import mpi4py.MPI
mpi4py.MPI.COMM_WORLD.Barrier()
start_t = time.time()
rss = []

sim.run(SIM_TIME)

rss.append(psutil.Process().memory_info().rss / 1024**2 - mem_ini)

rss = np.mean(rss)
avg_rss = mpi4py.MPI.COMM_WORLD.allreduce(rss, mpi4py.MPI.SUM) / MPI.nhosts
avg_mem_ini = mpi4py.MPI.COMM_WORLD.allreduce(mem_ini, mpi4py.MPI.SUM) / MPI.nhosts

mpi4py.MPI.COMM_WORLD.Barrier()
if MPI.rank == 0:
    if USE_STEPS_4:
        print("steps4 nhost %i time cost: %f" % (MPI.nhosts, time.time() - start_t))
        print("steps4 nhost %i avg rss (MB): %f" % (MPI.nhosts, avg_rss))
        print("steps4 nhost %i avg mem ini (MB): %f" % (MPI.nhosts, avg_mem_ini))
    else:
        print("steps3 nhost %i time cost: %f" % (MPI.nhosts, time.time() - start_t))
        print("steps3 nhost %i avg rss (MB): %f" % (MPI.nhosts, avg_rss))
        print("steps3 nhost %i avg mem ini (MB): %f" % (MPI.nhosts, avg_mem_ini))
