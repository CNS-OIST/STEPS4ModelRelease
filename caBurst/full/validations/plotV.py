from matplotlib import pyplot as plt

seeds = range(2001, 2101) # Scale 1 (length 10µm), smooth is first 4µm
#seeds = range(3001, 3101) # Scale 0.5 (length 5µm), smooth is first 0.5µm
seeds = range(4001, 4101) # STEPS 3, Scale 0.5 (length 5µm), smooth is first 0.5µm
#seeds = range(1001, 1101) # STEPS 4, Scale 0.5 (length 5µm), smooth is first 0.5µm

STEPS_version = 3

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
plt.show()

