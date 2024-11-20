import numpy as np
import pickle

seeds = range(1, 1000)
vlen = 301 # The number of time points voltage is recorded at per run

v1s = np.zeros((len(seeds), vlen))
v2s = np.zeros((len(seeds), vlen))
v3s = np.zeros((len(seeds), vlen))
v4s = np.zeros((len(seeds), vlen))
ts = None

for STEPS_version in [3, 4]:
    
    ofile1 = open(f'v1_{STEPS_version}.txt', 'wb')
    ofile2 = open(f'v2_{STEPS_version}.txt', 'wb')
    ofile3 = open(f'v3_{STEPS_version}.txt', 'wb')
    ofile4 = open(f'v4_{STEPS_version}.txt', 'wb')

    sidx = 0
    for seed in seeds:
        nd=np.genfromtxt(f"raw_traces/STEPS{STEPS_version}/respyramid_{seed}_STEPS{STEPS_version}.txt", delimiter=' ')[1:]
        
        v1s[sidx] = nd[:,0]
        v2s[sidx] = nd[:,1]
        v3s[sidx] = nd[:,2]
        v4s[sidx] = nd[:,3]
        
        ts = nd[:,-1]
        
        sidx += 1

    v1mean = np.mean(v1s, axis=0)
    v1std = np.std(v1s, axis=0)
    pickle.dump(v1mean, ofile1)
    pickle.dump(v1std, ofile1)
    
    v2mean = np.mean(v2s, axis=0)
    v2std = np.std(v2s, axis=0)
    pickle.dump(v2mean, ofile2)
    pickle.dump(v2std, ofile2)
    
    v3mean = np.mean(v3s, axis=0)
    v3std = np.std(v3s, axis=0)
    pickle.dump(v3mean, ofile3)
    pickle.dump(v3std, ofile3)

    v4mean = np.mean(v4s, axis=0)
    v4std = np.std(v4s, axis=0)
    pickle.dump(v4mean, ofile4)
    pickle.dump(v4std, ofile4)
    
    ofile1.close()
    ofile2.close()
    ofile3.close()
    ofile4.close()
