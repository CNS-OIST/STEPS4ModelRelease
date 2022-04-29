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

class SimulationError(Exception):
    pass

def run(seed, mesh_path, steps_version):
    if steps_version not in [3, 4]:
        raise SimulationError(f"Steps number: {steps_version} is not 3 or 4")

    # Potassium conductance, Siemens/m^2
    K_G = 360
    # Sodium conductance, Siemens/m^2
    Na_G = 1200
    # Leak conductance, Siemens/m^2
    L_G = 0.25

    # Potassium reversal potential, V
    K_rev = -77e-3
    # Sodium reversal potential, V
    Na_rev = 50e-3
    # Leak reveral potential, V
    leak_rev = -65.0e-3

    # Potassium channel density
    K_ro = 18.0e12
    # Sodium channel density
    Na_ro = 60.0e12

    # Total leak conductance for ideal cylinder:
    surfarea_cyl = 1.0*math.pi*1000*1e-12
    L_G_tot = L_G*surfarea_cyl


    # A table of potassium density factors at -65mV, found in getpops. n0, n1, n2, n3, n4
    K_FACS = [0.216750577045, 0.40366011853, 0.281904943772, \
                0.0874997924409, 0.0101845682113 ]

    # A table of sodium density factors. m0h1, m1h1, m2h1, m3h1, m0h0, m1h0, m2h0, m3h0
    NA_FACS = [0.343079175644, 0.0575250437508, 0.00321512825945, 5.98988373918e-05, \
                0.506380603793, 0.0849062503811, 0.00474548939393, 8.84099403236e-05]

    # Ohm.m
    Ra = 1.0

    # # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # #

    # Sim end time (seconds)
    SIM_END = 0.25

    # The current injection in amps
    Iinj = 0.1e-9

    EF_DT = 1e-6
    SAVE_DT = 5e-6

    # # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # #

    # record potential at the two extremes along (z) axis
    POT_POS = [0.0, 1.0e-03]

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Mesh geometry
    mesh = DistMesh(mesh_path) if steps_version == 4 else TetMesh.LoadGmsh(mesh_path)

    with mesh:
        if steps_version == 4:
            __MESH__ = Compartment.Create()

            memb = Patch.Create(__MESH__, None, 'ssys')
            z_min = Patch.Create(__MESH__, None, 'ssys')
            z_max = Patch.Create(__MESH__, None, 'ssys')
        else:
            __MESH__ = Compartment.Create(mesh.tets)

            memb = Patch.Create(mesh.triGroups[(0, 'memb')], __MESH__, None, 'ssys')
            z_min = Patch.Create(mesh.triGroups[(0, 'z_min')], __MESH__, None, 'ssys')
            z_max = Patch.Create(mesh.triGroups[(0, 'z_max')], __MESH__, None, 'ssys')


        surfarea_mesh = memb.Area
        corr_fac_area = surfarea_mesh/surfarea_cyl

        vol_cyl = math.pi*0.5*0.5*1000*1e-18
        vol_mesh = __MESH__.Vol
        corr_fac_vol = vol_mesh/vol_cyl

        if steps_version == 4:
            membrane = Membrane.Create([memb], capacitance=0.01 / corr_fac_area)
            __MESH__.Conductivity = 1 / (Ra * corr_fac_vol)
        else:
            membrane = Membrane.Create([memb])

    
        # The tetrahedrons from which to record potential
        POT_TET = TetList(mesh.tets[0, 0, z] for z in POT_POS)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    print("Create Model")

    mdl = Model()
    r = ReactionManager()
    with mdl:
        ssys = SurfaceSystem.Create()

        # K channel
        n0, n1, n2, n3, n4 = SubUnitState.Create()
        VGKC = Channel.Create([n0, n1, n2, n3, n4])

        # Na channel
        h0, h1, m0, m1, m2, m3 = SubUnitState.Create()
        mNaSU, hNaSU = SubUnit.Create(
            [m0, m1, m2, m3],
            [h0, h1],
        )
        VGNaC = Channel.Create([mNaSU, hNaSU])

        # Leak
        leaksus = SubUnitState.Create()
        Leak = Channel.Create([leaksus])

        # Gating kinetics
        a_n = VDepRate(lambda V: 1e3*((0.01*(10-(V*1e3+65))/(math.exp((10-(V*1e3+65))/10)-1))))
        b_n = VDepRate(lambda V: 1e3*((0.125*math.exp(-(V*1e3+65)/80))))
        a_m = VDepRate(lambda V: 1e3*((0.1*(25-(V*1e3+65))/(math.exp((25-(V*1e3+65))/10)-1))))
        b_m = VDepRate(lambda V: 1e3*((4*math.exp(-(V*1e3+65)/18))))
        a_h = VDepRate(lambda V: 1e3*((0.07*math.exp(-(V*1e3+65)/20))))
        b_h = VDepRate(lambda V: 1e3*((1/(math.exp((30-(V*1e3+65))/10)+1))))

        with ssys:
            with VGKC[...]:
                n0.s <r[1]> n1.s <r[2]> n2.s <r[3]> n3.s <r[4]> n4.s
                r[1].K = 4*a_n,   b_n
                r[2].K = 3*a_n, 2*b_n
                r[3].K = 2*a_n, 3*b_n
                r[4].K =   a_n, 4*b_n

            with VGNaC[...]:
                h1.s <r[1]> h0.s
                r[1].K = a_h, b_h

                m0.s <r[1]> m1.s <r[2]> m2.s <r[3]> m3.s
                r[1].K = 3*a_m,   b_m
                r[2].K = 2*a_m, 2*b_m
                r[3].K =   a_m, 3*b_m

            VGKC_I = OhmicCurr.Create(VGKC[n4], K_G / K_ro, K_rev)
            VGNaC_I = OhmicCurr.Create(VGNaC[m3, h0], Na_G / Na_ro, Na_rev)


            # Set the single-channel conductance:
            g_leak_sc = L_G_tot/len(membrane.tris)
            OC_L = OhmicCurr.Create(Leak[leaksus], g_leak_sc, leak_rev)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Create the solver objects
    rng = RNG('r123', 512, seed)
    #rng = RNG('mt19937', 512, seed)

    if steps_version == 4:
        sim = Simulation('DistTetOpSplit', mdl, mesh, rng)#, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
    else:
        part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
        sim = Simulation('TetOpSplit', mdl, mesh, rng, MPI.EF_DV_PETSC, part)

    # Data saving
    rs = ResultSelector(sim)

    Vrs = rs.TETS(POT_TET).V

    sim.toSave(Vrs, dt=SAVE_DT)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    sim.newRun()

    sim.TRIS(membrane.tris).Leak[leaksus].Count = 1

    sim.memb.VGNaC[m0, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[0])
    sim.memb.VGNaC[m1, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[1])
    sim.memb.VGNaC[m2, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[2])
    sim.memb.VGNaC[m3, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[3])
    sim.memb.VGNaC[m0, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[4])
    sim.memb.VGNaC[m1, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[5])
    sim.memb.VGNaC[m2, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[6])
    sim.memb.VGNaC[m3, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[7])
    sim.memb.VGKC[n0].Count = (K_ro*surfarea_cyl*K_FACS[0])
    sim.memb.VGKC[n1].Count = (K_ro*surfarea_cyl*K_FACS[1])
    sim.memb.VGKC[n2].Count = (K_ro*surfarea_cyl*K_FACS[2])
    sim.memb.VGKC[n3].Count = (K_ro*surfarea_cyl*K_FACS[3])
    sim.memb.VGKC[n4].Count = (K_ro*surfarea_cyl*K_FACS[4])

    sim.membrane.Potential = -65e-3

    minzverts = list(set([v for t in z_min.tris for v in t.verts]))
    for v in minzverts:
        sim.solver.setVertIClamp(v.idx, Iinj/len(minzverts))

    if steps_version != 4:
        sim.membrane.Capac = 0.01 / corr_fac_area
        sim.membrane.VolRes = Ra * corr_fac_vol

    sim.EfieldDT = EF_DT

    sim.run(SIM_END)

    if MPI.rank == 0:
        folder_path = os.path.join("raw_traces", f"STEPS{steps_version}")
        os.makedirs(folder_path, exist_ok=True)
        df = pd.DataFrame({"t":Vrs.time[0], "V_z_min":Vrs.data[0,:,0], "V_z_max":Vrs.data[0,:,1]})
        df.to_csv(folder_path + f"/res{seed}_STEPS{steps_version}.txt", sep=" ", index=False)

if __name__ == "__main__":
    run(seed=int(sys.argv[1]), mesh_path=sys.argv[2], steps_version=int(sys.argv[3]))

