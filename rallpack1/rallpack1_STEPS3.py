# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack1 for TetOpSplit

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os.path as osp
from random import *
import time

import sys

import numpy as np
import pandas as pd
import steps
from steps.geom import UNKNOWN_TET
import steps.geom as sgeom
import steps.model as smodel
import steps.rng as srng
import steps.mpi
import steps.mpi.solver as ssolver
import steps.utilities.meshio as meshio
import steps.utilities.geom_decompose as gd

class Rallpack1Params:
    # axial resistivity Ω·m
    Ra = 1.0
    # membrane resistivity Ω·m²
    Rm = 4.0
    # membrane capacity F/m²
    Cm = 0.01
    # p.d. across membrane V
    Em = -0.065
    # Leak reveral potential, V
    leak_rev = -65.0e-3

    # Total leak conductance for ideal cylinder:
    surfarea_cyl = 1.0*np.pi*1000*1e-12
    L_G_tot = surfarea_cyl/Rm

    # # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # #
    # efield solver
    SIM_EFIELD = ssolver.EF_DV_PETSC
    # The simulation dt (seconds); 
    SIM_DT = 5.0e-5
    # EFIELD DT
    SIM_EFIELD_DT = 1e-6
    # Sim end time (seconds)
    SIM_END = 0.000025
    # The number of sim 'time points'; * SIM_DT = sim end time
    SIM_NTPNTS = int(SIM_END/SIM_DT)+1

    # The current injection in amps
    Iinj = 0.1e-9

    # # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # #
    # record potential at the two extremes along (z) axis
    POT_POS = np.array([0.0, 1.0e-03])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def build_geometry(mesh_path, scale=1.):
    mesh, nodeproxy, tetproxy, triproxy = meshio.importGmsh(mesh_path, scale)

    trigroups = triproxy.getGroups()
    cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))
    memb_tris = trigroups[(0, "memb")]
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
    memb.addSurfsys('ssys')
    membrane = sgeom.Memb(
        'membrane', mesh, [memb], opt_method=2, search_percent=100.0)
    minzverts = []
    for tri in trigroups[(0, "z_min")]:
        minzverts.extend(mesh.getTri(tri))
    minzverts = list(set(minzverts))
    
    return mesh, minzverts, memb_tris

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def build_model(params, mesh, memb_tris):
    mdl = smodel.Model()
    ssys = smodel.Surfsys('ssys', mdl)

    # Leak
    L = smodel.Chan('L', mdl)
    Leak = smodel.ChanState('Leak', mdl, L)

    # Set the single-channel conductance:
    g_leak_sc = params.L_G_tot/len(memb_tris)
    OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate=Leak,
                            erev=params.leak_rev, g=g_leak_sc)

    return mdl
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def run_sim(seed, meshfile, meshscale, params):
    mesh, minzverts, memb_tris = build_geometry(meshfile, meshscale)
    if steps.mpi.rank == 0:
        meshio.saveMesh("compart", mesh)
    mdl = build_model(params, mesh, memb_tris)
    rng = srng.create('r123', 512)
    rng.initialize(seed)
    tet_hosts = gd.linearPartition(mesh, [1, 1, steps.mpi.nhosts])
    tri_hosts = gd.partitionTris(mesh, tet_hosts, memb_tris)

    # Create the solver objects
    sim = ssolver.TetOpSplit(
        mdl, mesh, rng, params.SIM_EFIELD, tet_hosts, tri_hosts)

    surfarea_mesh = sim.getPatchArea('memb')
    surfarea_cyl = 1.0*np.pi*1000*1e-12
    corr_fac_area = surfarea_mesh/surfarea_cyl

    vol_cyl = np.pi*0.5*0.5*1000*1e-18
    vol_mesh = sim.getCompVol('cyto')
    corr_fac_vol = vol_mesh/vol_cyl

    POT_N = len(params.POT_POS)
    POT_TET = []
    for p in params.POT_POS:
        POT_TET.append(mesh.findTetByPoint([0.0, 0.0, p]))
    RES_POT = np.zeros((params.SIM_NTPNTS, POT_N))

    for t in memb_tris:
        sim.setTriCount(t, 'Leak', 1)

    sim.setMembPotential('membrane', params.Em)
    sim.setMembVolRes('membrane', params.Ra*corr_fac_vol)
    sim.setMembCapac('membrane', params.Cm/corr_fac_area)

    for v in minzverts:
        sim.setVertIClamp(v, params.Iinj/len(minzverts))
    sim.setEfieldDT(params.SIM_EFIELD_DT)
    for l in range(params.SIM_NTPNTS):
        sim.run(params.SIM_DT*l)
        if steps.mpi.rank == 0:
            print(f"total progress: {round(1000*l/params.SIM_NTPNTS)/10.0}%")

        for p in range(POT_N):
            RES_POT[l, p] = sim.getTetV(int(POT_TET[p]))
    return RES_POT, np.array(range(params.SIM_NTPNTS))*params.SIM_DT, mesh

def gen_comparison(record, x0data, x1000data):
    v_benchmark_x0 = []
    v_benchmark_x1000 = []
    # At 0um- the end of the mesh
    with open(x0data) as istr:
        # Read in mv and ms
        next(istr)
        next(istr)
        for line in istr:
            nums = line.split()
            v_benchmark_x0.append(float(nums[1]))

    # At 1000um- the end of the mesh
    with open(x1000data) as istr:
        next(istr)
        next(istr)
        for line in istr:
            nums = line.split()
            v_benchmark_x1000.append(float(nums[1]))

    v_benchmark_x0 = v_benchmark_x0[:len(record)]
    v_benchmark_x1000 = v_benchmark_x1000[:len(record)]

    return [v_benchmark_x0, record[:, 0], v_benchmark_x1000, record[:, 1]]

def plot_data(simdata):
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.plot(simdata[0], 'k-' ,label = 'Benchmark, 0um', linewidth=3)
    plt.plot(simdata[1], 'r--', label = 'STEPS3 TetOpSplit, 0um', linewidth=3)
    plt.legend(loc='best')
    plt.ylabel('Potential (mV)')
    plt.subplot(212)
    plt.plot(simdata[2], 'k-' ,label = 'Benchmark, 1000um', linewidth=3)
    plt.plot(simdata[3], 'r--', label = 'STEPS3 TetOpSplit, 1000um', linewidth=3)
    plt.legend(loc='best')
    plt.ylabel('Potential (mV)')
    plt.xlabel('Time (ms)')
    plt.show()


def record_results(res_pot, t, mesh):
    df = pd.DataFrame({"t":t, "V_z_min":res_pot[:,0], "V_z_max":res_pot[:,1]})
    df.to_csv(f'results/STEPS3/res{int(sys.argv[1])}_STEPS3_{mesh.countVertices()*3}DoFs.txt', sep=" ", index=False)

if __name__ == "__main__":
    SEED = int(sys.argv[1])
    MESHFILE = sys.argv[2]
    # MESHFMT = 'msh'
    MESHSCALE = 1.0
    #X0DATA = '../benchmark_NEURON/v0'
    #X1000DATA = '../benchmark_NEURON/vx'
    sim_record, t, mesh = run_sim(SEED, MESHFILE, MESHSCALE, Rallpack1Params())
    #comparison_data = gen_comparison(sim_record, X0DATA, X1000DATA)
    #if steps.mpi.rank == 0:
    #    plot_data(comparison_data)
    record_results(sim_record, t, mesh)
