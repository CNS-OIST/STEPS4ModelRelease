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

# * : number of cores
for filename in sorted(glob.glob('strong_scaling_*.out'), key=numericalSort):    
    text = open(filename).read()
    text = text.split('\n')

    steps3_timings = []
    steps4_timings = []

    steps3_mem = []
    steps4_mem = []
    
    for line in text:
        if ("steps3" in line) and ("time cost" in line):
            line = line.split(" ")
            steps3_timings.append(float(line[-1]))
        if ("steps4" in line) and ("time cost" in line):
            line = line.split(" ")
            steps4_timings.append(float(line[-1]))
        
        if ("steps3" in line) and ("avg rss (MB)" in line):
            line = line.split(" ")
            steps3_mem.append(float(line[-1]))
        if ("steps4" in line) and ("avg rss (MB)" in line):
            line = line.split(" ")
            steps4_mem.append(float(line[-1]))

    print(int(filename.split("_")[-1].split(".")[0]), \
          np.mean(steps3_timings), np.std(steps3_timings), np.mean(steps3_mem), np.std(steps3_mem), \
          np.mean(steps4_timings), np.std(steps4_timings), np.mean(steps4_mem), np.std(steps4_mem))
