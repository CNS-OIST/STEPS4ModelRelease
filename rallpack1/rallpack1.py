# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# # Tests adapted from rallpack 1
# # Original Rallpack 1 author: Iain Hepburn
# # Test suite author: Alessandro Cattabiani
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function, absolute_import
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

class SimulationError(Exception):
    pass

def run(seed, mesh_path, steps_version):
    if steps_version not in [3, 4]:
        raise SimulationError(f"Steps number: {steps_version} is not 3 or 4")
    
    """Run rallpack 1 simulation"""

    # # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # #

    # Sim end time (seconds)
    SIM_END = 0.25

    # The current injection in amps
    Iinj = 0.1e-9

    EF_DT = 1e-6
    SAVE_DT = 5e-6

    # # # # # # # # # # # # # # # # PARAMETERS # # # # # # # # # # # # # #

    # Leak conductance, Siemens/m^2
    L_G = 0.25

    # Leak reveral potential, V
    leak_rev = -65.0e-3

    # Total leak conductance for ideal cylinder:
    surfarea_cyl = 1.0 * math.pi * 1000 * 1e-12
    L_G_tot = L_G * surfarea_cyl

    # Ohm.m
    Ra = 1.0

    # # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # #

    # record potential at the two extremes along (z) axis
    POT_POS = [0.0, 1.0e-03]

    # Mesh geometry
    mesh = DistMesh(mesh_path) if steps_version == 4 else TetMesh.LoadGmsh(mesh_path)

    with mesh:
        if steps_version == 4:
            __MESH__ = Compartment.Create()

            memb = Patch.Create(__MESH__, None, "ssys")
            z_min = Patch.Create(__MESH__, None, "ssys")
            z_max = Patch.Create(__MESH__, None, "ssys")
        else:
            __MESH__ = Compartment.Create(mesh.tets)

            memb = Patch.Create(mesh.triGroups[(0, "memb")], __MESH__, None, "ssys")
            z_min = Patch.Create(mesh.triGroups[(0, "z_min")], __MESH__, None, "ssys")
            z_max = Patch.Create(mesh.triGroups[(0, "z_max")], __MESH__, None, "ssys")

        surfarea_mesh = memb.Area
        corr_fac_area = surfarea_mesh / surfarea_cyl

        vol_cyl = math.pi * 0.5 * 0.5 * 1000 * 1e-18
        vol_mesh = __MESH__.Vol
        corr_fac_vol = vol_mesh / vol_cyl

        if steps_version == 4:
            membrane = Membrane.Create([memb], capacitance=0.01 / corr_fac_area)
            __MESH__.Conductivity = 1 / (Ra * corr_fac_vol)
        else:
            membrane = Membrane.Create([memb])

        # The tetrahedrons from which to record potential
        POT_TET = TetList(mesh.tets[0, 0, z] for z in POT_POS)

    # Model # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    mdl = Model()
    r = ReactionManager()
    with mdl:
        ssys = SurfaceSystem.Create()

        # Leak
        leaksus = SubUnitState.Create()

        Leak = Channel.Create([leaksus])

        with ssys:
            # Set the single-channel conductance:
            g_leak_sc = L_G_tot / len(membrane.tris)
            OC_L = OhmicCurr.Create(Leak[leaksus], g_leak_sc, leak_rev)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Create the solver objects
    rng = RNG("r123", 512, seed)

    if steps_version == 4:
        sim = Simulation(
            "DistTetOpSplit",
            mdl,
            mesh,
            rng,
            searchMethod=NextEventSearchMethod.GIBSON_BRUCK
        )
    else:
        part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
        sim = Simulation("TetOpSplit", mdl, mesh, rng, MPI.EF_DV_PETSC, part)

    # Data saving
    rs = ResultSelector(sim)

    Vrs = rs.TETS(POT_TET).V

    sim.toSave(Vrs, dt=SAVE_DT)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    sim.newRun()

    sim.TRIS(membrane.tris).Leak[leaksus].Count = 1

    sim.membrane.Potential = -65e-3

    minzverts = list(set([v for t in z_min.tris for v in t.verts]))
    for v in minzverts:
        sim.solver.setVertIClamp(v.idx, Iinj / len(minzverts))

    if steps_version != 4:
        sim.membrane.Capac = 0.01 / corr_fac_area
        sim.membrane.VolRes = Ra * corr_fac_vol

    sim.EfieldDT = EF_DT

    progress_dt = SAVE_DT * 10
    nsteps = math.floor(SIM_END / progress_dt) + 1

    for i in range(nsteps):
        current_time = i * progress_dt
        sim.run(current_time)
        print(
            f"Progress: {round(1e4 * (i / nsteps)) / 1e2}%, current time/SIM_END: "
            f"{round(1e4*current_time)/1e4}/{SIM_END}"
        )

    """Record"""
    folder_path = os.path.join("raw_traces", f"STEPS{steps_version}")
    os.makedirs(folder_path, exist_ok=True)
    df = pd.DataFrame(
        {"t": Vrs.time[0], "V zmin": Vrs.data[0, :, 0], "V zmax": Vrs.data[0, :, 1]}
    )

    df.to_csv(
        folder_path + f"/res{seed}_STEPS{steps_version}.txt",
        sep=" ",
        index=False,
    )

    return df


if __name__ == "__main__":
    run(seed=int(sys.argv[1]), mesh_path=sys.argv[2], steps_version=int(sys.argv[3]))
