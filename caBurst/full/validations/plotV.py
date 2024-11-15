from matplotlib import pyplot as plt
import numpy as np
import pickle

STEPS_version = 4

ofile1 = open(f'v1_{STEPS_version}.txt', 'wb')
ofile2 = open(f'v2_{STEPS_version}.txt', 'wb')
ofile3 = open(f'v3_{STEPS_version}.txt', 'wb')
ofile4 = open(f'v4_{STEPS_version}.txt', 'wb')

seeds = range(2001, 2101) # Scale 1 (length 10µm), smooth is first 4µm
#seeds = range(3001, 3101) # Scale 0.5 (length 5µm), smooth is first 0.5µm
seeds = range(4001, 4101) # STEPS 3, Scale 0.5 (length 5µm), smooth is first 0.5µm
#seeds = range(1001, 1101) # STEPS 4, Scale 0.5 (length 5µm), smooth is first 0.5µm
seeds = range(1, 1000)


v1s = [None]*len(seeds)
v2s = [None]*len(seeds)
v3s = [None]*len(seeds)
v4s = [None]*len(seeds)
ts = None
counter=0

v1s_5test = [None]*10
v2s_5test = [None]*10
v3s_5test = [None]*10
v4s_5test = [None]*10

counter5=0

for seed in seeds:
    ipath = f"raw_traces/STEPS{STEPS_version}/respyramid_{seed}_STEPS{STEPS_version}.txt"
    ifile = open(ipath)

    lines = ifile.readlines()

    v1=[]
    v2=[]
    v3=[]
    v4=[]
    t=[]

    for line in lines[1:]:
        line=line.split()
        v1.append(float(line[0]))
        v2.append(float(line[1]))
        v3.append(float(line[2]))
        v4.append(float(line[3]))
        t.append(float(line[-1]))
    
    v1s[counter] = v1
    v2s[counter] = v2
    v3s[counter] = v3
    v4s[counter] = v4
    ts=t

    counter+=1

    v1s_5test[counter5] = v1
    v2s_5test[counter5] = v2
    v3s_5test[counter5] = v3
    v4s_5test[counter5] = v4

    counter5+=1

    if not counter5%10:
        
        plt.subplot(221)
        v1mean = np.mean(v1s_5test, axis=0)
        plt.plot(t, v1mean, 'k-', linewidth=5)

        plt.subplot(222)
        v2mean = np.mean(v2s_5test, axis=0)
        plt.plot(t, v2mean, 'k-', linewidth=5)

        plt.subplot(223)
        v3mean = np.mean(v3s_5test, axis=0)
        plt.plot(t, v3mean, 'k-', linewidth=5)

        plt.subplot(224)
        v4mean = np.mean(v4s_5test, axis=0)
        plt.plot(t, v4mean, 'k-', linewidth=5)

        counter5 =0
            
    
    
    plt.subplot(221)
    plt.title("point1")
    plt.plot(t, v1, 'c-', lw=0.5, alpha=0.1)
    plt.ylim(-0.075, 0.02)
    plt.subplot(222)
    plt.title("point2")
    plt.plot(t, v2, 'm-', lw=0.5, alpha=0.1)
    plt.ylim(-0.075, 0.02)
    plt.subplot(223)
    plt.title("point3")
    plt.plot(t, v3, 'r-', lw=0.5, alpha=0.1)
    plt.ylim(-0.075, 0.02)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.subplot(224)
    plt.title("point4")
    plt.plot(t, v4, 'g-', lw=0.5, alpha=0.1)
    plt.ylim(-0.075, 0.02)
#plt.show()

#v1s=np.array(v1s)

plt.subplot(221)
v1mean = np.mean(v1s, axis=0)
v1std = np.std(v1s, axis=0)
plt.errorbar(t, v1mean, 2.0*v1std, alpha=0.5)
#plt.plot(t, v1mean, 'k-', linewidth=5)

plt.subplot(222)
v2mean = np.mean(v2s, axis=0)
v2std = np.std(v2s, axis=0)
plt.errorbar(t, v2mean, 2.0*v2std, alpha=0.5)
#plt.plot(t, v2mean, 'k-', linewidth=5)

plt.subplot(223)
v3mean = np.mean(v3s, axis=0)
v3std = np.std(v3s, axis=0)
plt.errorbar(t, v3mean, 2.0*v3std, alpha=0.5)
#plt.plot(t, v3mean, 'k-', linewidth=5)

plt.subplot(224)
v4mean = np.mean(v4s, axis=0)
v4std = np.std(v4s, axis=0)
plt.errorbar(t, v4mean, 2.0*v4std, alpha=0.5)
#plt.plot(t, v4mean, 'k-', linewidth=5)


plt.show()

pickle.dump(v1mean, ofile1)
pickle.dump(v1std, ofile1)

pickle.dump(v2mean, ofile2)
pickle.dump(v2std, ofile2)

pickle.dump(v3mean, ofile3)
pickle.dump(v3std, ofile3)

pickle.dump(v4mean, ofile4)
pickle.dump(v4std, ofile4)

ofile1.close()
ofile2.close()
ofile3.close()
ofile4.close()
