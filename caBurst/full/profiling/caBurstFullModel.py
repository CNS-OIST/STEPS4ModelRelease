########################################################################                                                                          
import steps.interface
                                                                                                                                                  
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

from extra.constants_withampa_yunliang import *

import numpy as np
import pandas as pd
import sys

import psutil
mem_ini = psutil.Process().memory_info().rss / 1024**2

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

BK_ro=BK_ro*1.5

SEED = int(sys.argv[1])
mesh_file = sys.argv[2]
USE_STEPS_4 = int(sys.argv[3]) > 0

######Glutamate transient#######
# Reference (Rudolph et al. 2011)
#Units (mM)
with open("../extra/Glut_Pulse_MVR.dat", "r") as f:
    Glut = list(map(float, f.read().split()))

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
r = ReactionManager()
with mdl:
    # Species
    Ca = Species.Create(valence=2)
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

    # Channels
    CaP_m0, CaP_m1, CaP_m2, CaP_m3 = SubUnitState.Create()
    CaPchan = Channel.Create([CaP_m0, CaP_m1, CaP_m2, CaP_m3])

    BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
    BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])

    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

    AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O = SubUnitState.Create()
    AMPA = Channel.Create([AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O])

    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])
    
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

    with ssys:
        # Ca Pump
        Pump.s + Ca.i <r[1]> CaPump.s >r[2]> Pump.s
        r[1].K = P_f_kcst, P_b_kcst
        r[2].K = P_k_kcst
    
        # CaP channel
        with CaPchan[...]:
            CaP_m0.s <r[1]> CaP_m1.s <r[2]> CaP_m2.s <r[3]> CaP_m3.s
            r[1].K = 3 * VDepRate(alpha_cap), 1 * VDepRate(beta_cap)
            r[2].K = 2 * VDepRate(alpha_cap), 2 * VDepRate(beta_cap)
            r[3].K = 1 * VDepRate(alpha_cap), 3 * VDepRate(beta_cap)
        OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, computeflux=True, virtual_oconc=Ca_oconc)
    
        # BK channel
        with BKchan[...]:
            (((BK_C0.s + Ca.i <r[1]> BK_C1.s)\
                       + Ca.i <r[2]> BK_C2.s)\
                       + Ca.i <r[3]> BK_C3.s)\
                       + Ca.i <r[4]> BK_C4.s
            r[1].K = c_01, c_10
            r[2].K = c_12, c_21
            r[3].K = c_23, c_32
            r[4].K = c_34, c_43

            (((BK_O0.s + Ca.i <r[1]> BK_O1.s)\
                       + Ca.i <r[2]> BK_O2.s)\
                       + Ca.i <r[3]> BK_O3.s)\
                       + Ca.i <r[4]> BK_O4.s
            r[1].K = o_01, o_10
            r[2].K = o_12, o_21
            r[3].K = o_23, o_32
            r[4].K = o_34, o_43
            
            BK_C0.s <r[1]> BK_O0.s
            BK_C1.s <r[2]> BK_O1.s
            BK_C2.s <r[3]> BK_O2.s
            BK_C3.s <r[4]> BK_O3.s
            BK_C4.s <r[5]> BK_O4.s
            r[1].K = VDepRate(f_0), VDepRate(b_0)
            r[2].K = VDepRate(f_1), VDepRate(b_1)
            r[3].K = VDepRate(f_2), VDepRate(b_2)
            r[4].K = VDepRate(f_3), VDepRate(b_3)
            r[5].K = VDepRate(f_4), VDepRate(b_4)
        OC_BK0 = OhmicCurr.Create(BKchan[BK_O0], BK_G, BK_rev)
        OC_BK1 = OhmicCurr.Create(BKchan[BK_O1], BK_G, BK_rev)
        OC_BK2 = OhmicCurr.Create(BKchan[BK_O2], BK_G, BK_rev)
        OC_BK3 = OhmicCurr.Create(BKchan[BK_O3], BK_G, BK_rev)
        OC_BK4 = OhmicCurr.Create(BKchan[BK_O4], BK_G, BK_rev)

        # SK channel
        with SKchan[...]:
            ((SK_C1.s + Ca.i <r[1]> SK_C2.s)\
                      + Ca.i <r[2]> SK_C3.s)\
                      + Ca.i <r[3]> SK_C4.s
            r[1].K = dirc2_t, invc1_t
            r[2].K = dirc3_t, invc2_t
            r[3].K = dirc4_t, invc3_t
            
            SK_C3.s <r[1]> SK_O1.s
            SK_C4.s <r[2]> SK_O2.s
            r[1].K = diro1_t, invo1_t
            r[2].K = diro2_t, invo2_t
        OC1_SK = OhmicCurr.Create(SKchan[SK_O1], SK_G, SK_rev)
        OC2_SK = OhmicCurr.Create(SKchan[SK_O2], SK_G, SK_rev)

        # AMPA channel
        with AMPA[...]:
            AMPA_C.s <r['AMPACC1']> AMPA_C1.s <r['AMPAC1C2']> AMPA_C2.s <r[3]> AMPA_O.s
            r['AMPACC1'].K = 0, ru1
            r['AMPAC1C2'].K = 0, ru2
            r[3].K = ro, rc

            AMPA_C1.s <r[1]> AMPA_D1.s
            AMPA_C2.s <r[2]> AMPA_D2.s
            r[1].K = rd, rr
            r[2].K = rd, rr
        OC_AMPAR1 = OhmicCurr.Create(AMPA[AMPA_O], AMPA_G, AMPA_rev)

        # Leak current channel
        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)
    
########### MESH & COMPARTMENTALIZATION #################

mesh = DistMesh(mesh_file, 1e-6) if USE_STEPS_4 else TetMesh.LoadGmsh(mesh_file, 1e-6)

with mesh:
    if USE_STEPS_4:
        __MESH__ = Compartment.Create(vsys, conductivity=1 / Ra)

        smooth = Patch(__MESH__, None, ssys, name='smooth.__BOUNDARY__')
        spiny = Patch(__MESH__, None, ssys, name='spiny.__BOUNDARY__')

        memb_spiny = Membrane.Create([spiny], capacitance=memb_capac_spiny)
        memb_smooth = Membrane.Create([smooth], capacitance=memb_capac_proximal)
    else:
        __MESH__ = Compartment.Create(mesh.tets, vsys)

        smoothTris = mesh.tetGroups[(0, 'smooth')].surface & mesh.surface
        spinyTris = mesh.tetGroups[(0, 'spiny')].surface & mesh.surface
        smooth = Patch(smoothTris, __MESH__, None, ssys, name='smooth.__BOUNDARY__')
        spiny = Patch(spinyTris, __MESH__, None, ssys, name='spiny.__BOUNDARY__')

        memb = Membrane.Create([spiny, smooth])

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, SEED)

if USE_STEPS_4:
    sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
else:
    part = GmshPartition(mesh)
    sim = Simulation('TetOpSplit', mdl, mesh, rng, MPI.EF_DV_PETSC, part)

# set temperature for ghk reactions
sim.Temp = 273.15+TEMPERATURE

if MPI.rank == 0:
    print("solver created.", flush=True)

rs = ResultSelector(sim)

Pots = rs.MAX(rs.TRIS(smooth.tris).V) << \
       rs.MIN(rs.TRIS(smooth.tris).V) << \
       rs.MAX(rs.TRIS(spiny.tris).V) << \
       rs.MIN(rs.TRIS(spiny.tris).V)

sim.toSave(Pots, dt=TIMECONVERTER * 10)

sim.newRun()

smooth_area = smooth.Area
if MPI.rank == 0:
    print("smooth area: ", smooth_area)
spiny_area = spiny.Area
if MPI.rank == 0:
    print("spiny area: ", spiny_area)

#Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
#pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
pumpnbs = 6.022141e12*smooth_area

sim.LIST(smooth.name).Pump.Count = round(pumpnbs)
sim.LIST(smooth.name).CaPump.Count = 0

sim.LIST(smooth.name).CaPchan[CaP_m0].Count = round(CaP_ro*smooth_area*CaP_m0_p)
sim.LIST(smooth.name).CaPchan[CaP_m1].Count = round(CaP_ro*smooth_area*CaP_m1_p)
sim.LIST(smooth.name).CaPchan[CaP_m2].Count = round(CaP_ro*smooth_area*CaP_m2_p)
sim.LIST(smooth.name).CaPchan[CaP_m3].Count = round(CaP_ro*smooth_area*CaP_m3_p)

sim.LIST(smooth.name).BKchan[BK_C0].Count = round(BK_ro*smooth_area*BK_C0_p)
sim.LIST(smooth.name).BKchan[BK_C1].Count = round(BK_ro*smooth_area*BK_C1_p)
sim.LIST(smooth.name).BKchan[BK_C2].Count = round(BK_ro*smooth_area*BK_C2_p)
sim.LIST(smooth.name).BKchan[BK_C3].Count = round(BK_ro*smooth_area*BK_C3_p)
sim.LIST(smooth.name).BKchan[BK_C4].Count = round(BK_ro*smooth_area*BK_C4_p)

sim.LIST(smooth.name).BKchan[BK_O0].Count = round(BK_ro*smooth_area*BK_O0_p)
sim.LIST(smooth.name).BKchan[BK_O1].Count = round(BK_ro*smooth_area*BK_O1_p)
sim.LIST(smooth.name).BKchan[BK_O2].Count = round(BK_ro*smooth_area*BK_O2_p)
sim.LIST(smooth.name).BKchan[BK_O3].Count = round(BK_ro*smooth_area*BK_O3_p)
sim.LIST(smooth.name).BKchan[BK_O4].Count = round(BK_ro*smooth_area*BK_O4_p)

sim.LIST(smooth.name).SKchan[SK_C1].Count = round(SK_ro*smooth_area*SK_C1_p)
sim.LIST(smooth.name).SKchan[SK_C2].Count = round(SK_ro*smooth_area*SK_C2_p)
sim.LIST(smooth.name).SKchan[SK_C3].Count = round(SK_ro*smooth_area*SK_C3_p)
sim.LIST(smooth.name).SKchan[SK_C4].Count = round(SK_ro*smooth_area*SK_C4_p)

sim.LIST(smooth.name).SKchan[SK_O1].Count = round(SK_ro*smooth_area*SK_O1_p)
sim.LIST(smooth.name).SKchan[SK_O2].Count = round(SK_ro*smooth_area*SK_O2_p)

sim.LIST(smooth.name).AMPA[AMPA_C].Count = round(AMPA_receptors)
sim.LIST(smooth.name).AMPA[AMPA_C1].Count = 0
sim.LIST(smooth.name).AMPA[AMPA_C2].Count = 0
sim.LIST(smooth.name).AMPA[AMPA_O].Count = 0
sim.LIST(smooth.name).AMPA[AMPA_D1].Count = 0
sim.LIST(smooth.name).AMPA[AMPA_D2].Count = 0

sim.LIST(smooth.name).L[Leak].Count = round(L_ro_proximal * smooth_area)

#Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
#pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
pumpnbs = 6.022141e12*spiny_area

sim.LIST(spiny.name).Pump.Count = round(pumpnbs)
sim.LIST(spiny.name).CaPump.Count = 0

sim.LIST(spiny.name).CaPchan[CaP_m0].Count = round(CaP_ro*spiny_area*CaP_m0_p)
sim.LIST(spiny.name).CaPchan[CaP_m1].Count = round(CaP_ro*spiny_area*CaP_m1_p)
sim.LIST(spiny.name).CaPchan[CaP_m2].Count = round(CaP_ro*spiny_area*CaP_m2_p)
sim.LIST(spiny.name).CaPchan[CaP_m3].Count = round(CaP_ro*spiny_area*CaP_m3_p)

sim.LIST(spiny.name).BKchan[BK_C0].Count = round(BK_ro*spiny_area*BK_C0_p)
sim.LIST(spiny.name).BKchan[BK_C1].Count = round(BK_ro*spiny_area*BK_C1_p)
sim.LIST(spiny.name).BKchan[BK_C2].Count = round(BK_ro*spiny_area*BK_C2_p)
sim.LIST(spiny.name).BKchan[BK_C3].Count = round(BK_ro*spiny_area*BK_C3_p)
sim.LIST(spiny.name).BKchan[BK_C4].Count = round(BK_ro*spiny_area*BK_C4_p)

sim.LIST(spiny.name).BKchan[BK_O0].Count = round(BK_ro*spiny_area*BK_O0_p)
sim.LIST(spiny.name).BKchan[BK_O1].Count = round(BK_ro*spiny_area*BK_O1_p)
sim.LIST(spiny.name).BKchan[BK_O2].Count = round(BK_ro*spiny_area*BK_O2_p)
sim.LIST(spiny.name).BKchan[BK_O3].Count = round(BK_ro*spiny_area*BK_O3_p)
sim.LIST(spiny.name).BKchan[BK_O4].Count = round(BK_ro*spiny_area*BK_O4_p)

sim.LIST(spiny.name).SKchan[SK_C1].Count = round(SK_ro*spiny_area*SK_C1_p)
sim.LIST(spiny.name).SKchan[SK_C2].Count = round(SK_ro*spiny_area*SK_C2_p)
sim.LIST(spiny.name).SKchan[SK_C3].Count = round(SK_ro*spiny_area*SK_C3_p)
sim.LIST(spiny.name).SKchan[SK_C4].Count = round(SK_ro*spiny_area*SK_C4_p)

sim.LIST(spiny.name).SKchan[SK_O1].Count = round(SK_ro*spiny_area*SK_O1_p)
sim.LIST(spiny.name).SKchan[SK_O2].Count = round(SK_ro*spiny_area*SK_O2_p)

sim.LIST(spiny.name).AMPA[AMPA_C].Count = 0
sim.LIST(spiny.name).AMPA[AMPA_C1].Count = 0
sim.LIST(spiny.name).AMPA[AMPA_C2].Count = 0
sim.LIST(spiny.name).AMPA[AMPA_O].Count = 0
sim.LIST(spiny.name).AMPA[AMPA_D1].Count = 0
sim.LIST(spiny.name).AMPA[AMPA_D2].Count = 0

sim.LIST(spiny.name).L[Leak].Count = round(L_ro_spiny * spiny_area)

sim.__MESH__.Ca.Conc = Ca_iconc
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

if not USE_STEPS_4:
    sim.ALL(Membrane).VolRes = Ra
    sim.TRIS(smooth.tris).Capac = memb_capac_proximal
    sim.TRIS(spiny.tris).Capac = memb_capac_spiny

sim.EfieldDT = EF_DT

sim.ALL(Membrane).Potential = init_pot


if MPI.rank == 0:
    print("simulation start")

NTIMEPOINTS =  21 #3001 | found in extra/constants_withampa_yunliang.py

import time
import mpi4py.MPI
mpi4py.MPI.COMM_WORLD.Barrier()
start_t = time.time()
rss = []
for l in range(NTIMEPOINTS):
    t_now = time.time()
    if MPI.rank == 0:
        print("Tpnt: ", l, "/", NTIMEPOINTS)
        print("Sim Time: ", 1.0e3 * TIMECONVERTER * l)

    sim.LIST(smooth.name).AMPACC1['fwd'].K = 1.0e-3 * rb * Glut[l + 2495]
    sim.LIST(smooth.name).AMPAC1C2['fwd'].K = 1.0e-3 * rb * Glut[l + 2495]

    sim.run(TIMECONVERTER * l)

    rss.append(psutil.Process().memory_info().rss / 1024**2 - mem_ini)
    if MPI.rank == 0:
        print("*-*-*->", "sim time per iT (sec) : ", time.time() - t_now, flush=True)
        print("*-*-*->", "iT : Done", flush=True) 

rss = np.mean(rss)
avg_rss = mpi4py.MPI.COMM_WORLD.allreduce(rss, mpi4py.MPI.SUM) / MPI.nhosts
avg_mem_ini = mpi4py.MPI.COMM_WORLD.allreduce(mem_ini, mpi4py.MPI.SUM) / MPI.nhosts

mpi4py.MPI.COMM_WORLD.Barrier()
if MPI.rank == 0:
    print("end simulation")
    print("start recording")
    labels = ["smooth_max_V", "smooth_min_V", "spiny_max_V", "spiny_min_V"]
    dct = {name: Pots.data[0,:,i] for i, name in enumerate(labels)}
    dct["t"] = Pots.time[0]
    pd.DataFrame(dct).to_csv(f'results/{MPI.nhosts}_res{SEED}_STEPS{4 if USE_STEPS_4 else 3}.txt', sep=" ", index=False)

if MPI.rank == 0:
    if USE_STEPS_4:
        print("steps4 nhost %i time cost: %f" % (MPI.nhosts, time.time() - start_t))
        print("steps4 nhost %i avg rss (MB): %f" % (MPI.nhosts, avg_rss))
        print("steps4 nhost %i avg mem ini (MB): %f" % (MPI.nhosts, avg_mem_ini))
    else:
        print("steps3 nhost %i time cost: %f" % (MPI.nhosts, time.time() - start_t))
        print("steps3 nhost %i avg rss (MB): %f" % (MPI.nhosts, avg_rss))
        print("steps3 nhost %i avg mem ini (MB): %f" % (MPI.nhosts, avg_mem_ini))
