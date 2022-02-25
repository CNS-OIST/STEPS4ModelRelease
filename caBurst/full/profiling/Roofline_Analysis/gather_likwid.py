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

for filename in sorted(glob.glob('likwid_*.out'), key=numericalSort):
    print(filename)

print(80*"*")

lst = 0
for filename in sorted(glob.glob('likwid_*.out'), key=numericalSort):
    text = open(filename).read()

    text = text.split('\n')

    regions = []
    DP = {}
    BW = {}
    for line in text:
        if "Region:" in line:
            regions.append(line)
        if "AVX DP [MFLOP/s] STAT" in line:
            continue
        if "DP [MFLOP/s] STAT" in line:
            line = line.split("|")
            if (regions[-1] == "Region: SSAOperator::run()") or \
               (regions[-1] == "Region: DiffusionOperator::operator()") or \
               (regions[-1] == "Region: EFieldOperator::evolve()"):
               DP[regions[-1]] = float(line[2])
        if "Memory bandwidth [MBytes/s] STAT" in line:
            line = line.split("|")
            if (regions[-1] == "Region: SSAOperator::run()") or \
               (regions[-1] == "Region: DiffusionOperator::operator()") or \
               (regions[-1] == "Region: EFieldOperator::evolve()"):
               BW[regions[-1]] = float(line[2])

    TASKS = int(filename.split(".")[0].split("_")[-1])
    _1 = DP["Region: SSAOperator::run()"]
    _2 = BW["Region: SSAOperator::run()"]
    _3 = DP["Region: DiffusionOperator::operator()"]
    _4 = BW["Region: DiffusionOperator::operator()"]
    _5 = DP["Region: EFieldOperator::evolve()"]
    _6 = BW["Region: EFieldOperator::evolve()"]
    print(f"node_spread[{lst}] = [{TASKS}, [[{_1}, {_2}], [{_3}, {_4}], [{_5}, {_6}]]]")
    lst += 1

print(80*"*")
