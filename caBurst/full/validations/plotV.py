from pylab import *

seeds = range(2001, 2101) # Scale 1, smooth is first 4e-6
seeds = range(3001, 3101) # Scale 0.5, smooth is first 1e-6
for seed in seeds:
    ipath = f"raw_traces/STEPS3/respyramid_{seed}_STEPS3.txt"
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
    
    subplot(221)
    title("point1")
    plot(t, v1, 'c-', lw=0.5, alpha=0.1)
    ylim(-0.075, 0.02)
    subplot(222)
    title("point2")
    plot(t, v2, 'm-', lw=0.5, alpha=0.1)
    ylim(-0.075, 0.02)
    subplot(223)
    title("point3")
    plot(t, v3, 'r-', lw=0.5, alpha=0.1)
    ylim(-0.075, 0.02)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    subplot(224)
    title("point4")
    plot(t, v4, 'g-', lw=0.5, alpha=0.1)
    ylim(-0.075, 0.02)
show()

