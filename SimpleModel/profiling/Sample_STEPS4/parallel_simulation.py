#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *
import time
import os
import sys

import numpy as np
import psutil
mem_ini = psutil.Process().memory_info().rss / 1024**2

if len(sys.argv) == 3:
    MOLECULE_RATIO = float(sys.argv[1])
    MESHFILE = sys.argv[2]
else:
    MOLECULE_RATIO = 1.0
    MESHFILE = '10x10x100_3363tets.msh'

# The initial molecule counts
N0A = 1000
N0B = 2000
N0C = 3000
N0D = 4000
N0E = 5000
N0F = 6000
N0G = 7000
N0H = 8000
N0I = 9000
N0J = 10000

ENDTIME = 0.5
DT = 0.001

########################################################################
# Biochemical Model
    
mdl = Model()
r = ReactionManager()
with mdl:
    # The chemical species
    SA, SB, SC, SD, SE, SF, SG, SH, SI, SJ = Species.Create()

    vsys = VolumeSystem.Create()

    with vsys:

        SA + SB <r[1]> SC
        SC + SD <r[2]> SE
        SF + SG <r[3]> SH
        SH + SI <r[4]> SJ
        r[1].K = 1000e6, 100
        r[2].K = 100e6, 10
        r[3].K = 10e6, 1
        r[4].K = 1e6, 1

        # The diffusion rules
        D1 = Diffusion.Create(SA, 100e-12)
        D2 = Diffusion.Create(SB, 90e-12)
        D3 = Diffusion.Create(SC, 80e-12)
        D4 = Diffusion.Create(SD, 70e-12)
        D5 = Diffusion.Create(SE, 60e-12)
        D6 = Diffusion.Create(SF, 50e-12)
        D7 = Diffusion.Create(SG, 40e-12)
        D8 = Diffusion.Create(SH, 30e-12)
        D9 = Diffusion.Create(SI, 20e-12)
        D10 = Diffusion.Create(SJ, 10e-12)

########################################################################
# Geometry

mesh = DistMesh(os.path.join('../mesh/', MESHFILE), scale=1e-6)
with mesh:
    comp1 = Compartment.Create(vsys, tetLst=mesh.tets)

########################################################################

rng = RNG('mt19937', 512, int(time.time()%4294967295))
sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)

rs = ResultSelector(sim)

counts = rs.comp1.ALL(Species).Count
counts.toFile('results.dat')

sim.toSave(counts, dt=DT)

########################################################################

sim.newRun()

# Set initial conditions
sim.comp1.SA.Count = N0A * MOLECULE_RATIO
sim.comp1.SB.Count = N0B * MOLECULE_RATIO
sim.comp1.SC.Count = N0C * MOLECULE_RATIO
sim.comp1.SD.Count = N0D * MOLECULE_RATIO
sim.comp1.SE.Count = N0E * MOLECULE_RATIO
sim.comp1.SF.Count = N0F * MOLECULE_RATIO
sim.comp1.SG.Count = N0G * MOLECULE_RATIO
sim.comp1.SH.Count = N0H * MOLECULE_RATIO
sim.comp1.SI.Count = N0I * MOLECULE_RATIO
sim.comp1.SJ.Count = N0J * MOLECULE_RATIO

import time
import mpi4py.MPI
mpi4py.MPI.COMM_WORLD.Barrier()
start_t = time.time()
rss = []
sim.run(ENDTIME)
rss.append(psutil.Process().memory_info().rss / 1024**2 - mem_ini)

rss = np.mean(rss)
avg_rss = mpi4py.MPI.COMM_WORLD.allreduce(rss, mpi4py.MPI.SUM) / MPI.nhosts
avg_mem_ini = mpi4py.MPI.COMM_WORLD.allreduce(mem_ini, mpi4py.MPI.SUM) / MPI.nhosts

mpi4py.MPI.COMM_WORLD.Barrier()
if MPI.rank == 0:
    print("steps4 nhost %i time cost: %f" % (MPI.nhosts, time.time() - start_t))
    print("steps4 nhost %i avg rss (MB): %f" % (MPI.nhosts, avg_rss))
    print("steps4 nhost %i avg mem ini (MB): %f" % (MPI.nhosts, avg_mem_ini))
