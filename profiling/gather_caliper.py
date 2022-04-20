import os
import glob
import re
import numpy as np

# For numerical ordering of files in a folder
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

# *.out : number of cores
for filename in sorted(glob.glob('*caliper_*.out'), key=numericalSort):    
    text = open(filename).read()
    text = text.split('\n')

    timings = [] # Wall Clock of STEPS4
    EField = []
    Diffusion = []
    SSA = []
    
    for line in text:
        if "time cost" in line:
            line = line.split(" ")
            timings.append(float(line[-1]))
        if "EFieldOperator::evolve()" in line:
            line = line.split(" ")
            line = [el for el in line if el !=""]
            EField.append(float(line[1])) # min time
            EField.append(float(line[2])) # max time
            EField.append(float(line[3])) # avg time
        if "DiffusionOperator::operator()" in line:
            line = line.split(" ")
            line = [el for el in line if el !=""]
            Diffusion.append(float(line[1]))
            Diffusion.append(float(line[2]))
            Diffusion.append(float(line[3]))
        if "SSAOperator::run()" in line:
            line = line.split(" ")
            line = [el for el in line if el !=""]
            SSA.append(float(line[1]))
            SSA.append(float(line[2]))
            SSA.append(float(line[3]))
    
    timings = np.mean(timings)
    EField = np.array([np.mean(EField[0::3]), np.mean(EField[1::3]), np.mean(EField[2::3])]) if len(EField) != 0 else [0.,0.,0.]
    Diffusion = np.array([np.mean(Diffusion[0::3]), np.mean(Diffusion[1::3]), np.mean(Diffusion[2::3])]) if len(Diffusion) != 0 else [0.,0.,0.]
    SSA = np.array([np.mean(SSA[0::3]), np.mean(SSA[1::3]), np.mean(SSA[2::3])]) if len(SSA) != 0 else [0.,0.,0.]

    PER_EField = (EField/timings)*100
    PER_Diffusion = (Diffusion/timings)*100
    PER_SSA = (SSA/timings)*100

    PER_OTHER = np.array([100.,100.,100.]) - PER_EField - PER_Diffusion - PER_SSA

    print(int(filename.split("_")[-1].split(".")[0]), \
          *PER_EField, *PER_Diffusion, *PER_SSA, *PER_OTHER, \
          *EField, *Diffusion, *SSA)
