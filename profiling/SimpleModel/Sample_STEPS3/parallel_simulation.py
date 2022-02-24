#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps
import steps.mpi
import steps.model as smod
import steps.utilities.geom_decompose as gd
import steps.mpi.solver as ssolver
import steps.utilities.meshio as meshio
import steps.geom as sgeom
import steps.rng as srng
import time
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

def gen_model():
    
    mdl = smod.Model()
    
    # The chemical species
    A = smod.Spec('A', mdl)
    B = smod.Spec('B', mdl)
    C = smod.Spec('C', mdl)
    D = smod.Spec('D', mdl)
    E = smod.Spec('E', mdl)
    F = smod.Spec('F', mdl)
    G = smod.Spec('G', mdl)
    H = smod.Spec('H', mdl)
    I = smod.Spec('I', mdl)
    J = smod.Spec('J', mdl)

    volsys = smod.Volsys('vsys',mdl)


    R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C],  kcst = 1000.0e6)
    R2 = smod.Reac('R2', volsys, lhs = [C],  rhs = [A,B], kcst = 100)
    R3 = smod.Reac('R3', volsys, lhs = [C, D], rhs = [E], kcst = 100e6)
    R4 = smod.Reac('R4', volsys, lhs = [E], rhs = [C,D], kcst = 10)

    R5 = smod.Reac('R5', volsys, lhs = [F, G], rhs = [H], kcst = 10e6)
    R6 = smod.Reac('R6', volsys, lhs = [H], rhs = [F,G], kcst = 1)
    R7 = smod.Reac('R7', volsys, lhs = [H, I], rhs = [J],  kcst = 1e6)
    R8 = smod.Reac('R8', volsys, lhs = [J],  rhs = [H,I], kcst = 1)


    # The diffusion rules
    D1 = smod.Diff('D1', volsys, A,  100e-12)
    D2 = smod.Diff('D2', volsys, B,  90e-12)
    D3 = smod.Diff('D3', volsys, C, 80e-12)
    D4 = smod.Diff('D4', volsys, D, 70e-12)
    D5 = smod.Diff('D5', volsys, E, 60e-12)
    D6 = smod.Diff('D6', volsys, F,  50e-12)
    D7 = smod.Diff('D7', volsys, G,  40e-12)
    D8 = smod.Diff('D8', volsys, H,  30e-12)
    D9 = smod.Diff('D9', volsys, I,  20e-12)
    D10 = smod.Diff('D10', volsys, J, 10e-12)
    
    return mdl

########################################################################
# Geometry
def gen_geom():
    mesh = meshio.importGmsh("../mesh/" + MESHFILE, 1e-6)[0]
    ntets = mesh.countTets()
    comp = sgeom.TmComp('comp', mesh, range(ntets))
    comp.addVolsys('vsys')
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

########################################################################
rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295)) # The max unsigned long

# Partitioning
tet_hosts = gd.linearPartition(g, [1,1,steps.mpi.nhosts])

sim = ssolver.TetOpSplit(m, g, rng, False, tet_hosts)

########################################################################

sim.reset()

# Set initial conditions
sim.setCompCount('comp', 'A', N0A * MOLECULE_RATIO)
sim.setCompCount('comp', 'B', N0B * MOLECULE_RATIO)
sim.setCompCount('comp', 'C', N0C * MOLECULE_RATIO)
sim.setCompCount('comp', 'D', N0D * MOLECULE_RATIO)
sim.setCompCount('comp', 'E', N0E * MOLECULE_RATIO)
sim.setCompCount('comp', 'F', N0F * MOLECULE_RATIO)
sim.setCompCount('comp', 'G', N0G * MOLECULE_RATIO)
sim.setCompCount('comp', 'H', N0H * MOLECULE_RATIO)
sim.setCompCount('comp', 'I', N0I * MOLECULE_RATIO)
sim.setCompCount('comp', 'J', N0J * MOLECULE_RATIO)

specNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

data = []

import time
import mpi4py.MPI
mpi4py.MPI.COMM_WORLD.Barrier()
start_t = time.time()
rss = []
for i in range(int(ENDTIME // DT) + 1):
    sim.run(i * DT)

    counts = [sim.getCompCount('comp', spec) for spec in specNames] 
    data.append([i * DT,] + counts)
    
    rss.append(psutil.Process().memory_info().rss / 1024**2 - mem_ini)

rss = np.mean(rss)
avg_rss = mpi4py.MPI.COMM_WORLD.allreduce(rss, mpi4py.MPI.SUM) / steps.mpi.nhosts
avg_mem_ini = mpi4py.MPI.COMM_WORLD.allreduce(mem_ini, mpi4py.MPI.SUM) / steps.mpi.nhosts

mpi4py.MPI.COMM_WORLD.Barrier()
if steps.mpi.rank == 0:
    with open('results.csv', 'w') as f:
        f.write(','.join(['t'] + specNames) + '\n')
        for row in data:
            f.write(','.join(map(str, row)) + '\n')
    
    print("steps3 nhost %i time cost: %f" % (steps.mpi.nhosts, time.time() - start_t))
    print("steps3 nhost %i avg rss (MB): %f" % (steps.mpi.nhosts, avg_rss))
    print("steps3 nhost %i avg mem ini (MB): %f" % (steps.mpi.nhosts, avg_mem_ini))
