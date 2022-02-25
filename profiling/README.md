Here, we describe how to perform the strong scaling & Roofline analysis found in STEPS 4 Release paper. To simply re-generate the graphs of the paper just run the `STEPS_PerfGraphs.ipynb` jupyter notebook (all data baked in).

## How to perform a Strong Scaling analysis

There are 3 main case studies, namely SimpleModel, caBurstBackground, and caBurstFull, for which we are interested in performing a strong scaling analysis. To perform a strong scaling analysis, the problem size remains fixed while the number of cores is steadily increased. The steps to do so are presented per case study below:

1. SimpleModel (see corresponding folder `../SimpleModel/profiling`):
    a. Submit multiple jobs in the queue of a supercomputer (in our case is BB5 of the Blue Brain Project) with the help of `strong_scaling.batch` script:
    ```
    sbatch --job-name=strong_scaling_2 --nodes=1 --ntasks-per-node=2 strong_scaling.batch
    ...
    sbatch --job-name=strong_scaling_2048 --nodes=64 --ntasks-per-node=32 strong_scaling.batch
    ...
    ```
    b. Once the jobs are finished, run the `strong_scaling.py` python script (found in this folder) to extract the needed data, and use it in the jupyter notebook. The corresponding cell in the `STEPS_PerfGraphs.ipynb` jupyter notebook generates the strong scaling graphs.

2. caBurstBackground (see corresponding folder `../caBurst/background/profiling`):
    a. Exact same steps as in **SimpleModel**

3. caBurstFull (see corresponding folder `../caBurst/full/profiling`):
    a. Exact same steps as in **SimpleModel**

### Strong Scaling analysis : Caliper Instrumentation

For every case study (namely SimpleModel, caBurstBackground, caBurstFull), there is a job submission script that uses Caliper for further instrumentation of the code, `CaseStudy/profiling/caliper.sbatch` (just for STEPS4). One could repeat the strong scaling analysis as described above, but using this script instead. Once all jobs per case study are finished, then just use the `gather_caliper.py` python script (found in this folder) to extract the data. This data per case study should be used in the `STEPS_PerfGraphs.ipynb` jupyter notebook, which generates the Caliper-based graphs. Keep in mind that STEPS must be compiled with Caliper support `ON` (see CMake variable `USE_CALIPER_PROFILING`).

## How to perform a Roofline analysis

We have chosen to perform the Roofline analysis for the caBurstFull case study given that it combines all the main computational kernels, i.e., SSA operator, Diffusion operator, and Efield operator. In the `../caBurst/full/profiling/Roofline_Analysis/` folder, one can find the job submission script (always single node analysis), and the python scipt to run the caBurstFull model:

* Submit multiple jobs in the queue of a supercomputer (in our case is BB5 of the Blue Brain Project) with the help of `likwid.batch` script:
    ```
    export NTASKS=2 
    sbatch --job-name=likwid_2 likwid.batch
    ...
    export NTASKS=32 
    sbatch --job-name=likwid_32 likwid.batch
    ...
    ```

* Once the jobs are finished, one can extract the data with the `gather_likwid.py` python script. The data should be used in the corresponding cell (last one) of the `STEPS_PerfGraphs.ipynb` jupyter notebook to generate the Roofline graph.

Few important things to highlight:
1. STEPS must be compiled with LIKWID support `ON` (see CMake variable `USE_LIKWID_PROFILING`).

2. The `caBurstFullModel.py` script found in `Roofline_Analysis/` folder initializes correctly the LIKWID environment. This is done through LIKWID's python bindings (see, `import pylikwid, etc.`).

3. It is very crucial to use `export STEPS_INSTRUMENTOR_MPI_BARRIER="before;after"` to make sure that we profile with LIKWID the annotated regions correctly. This happens in `likwid.batch` job submission script. The MPI barrier guarantees that all the processes traverse at the same time the annotated computational kernels, and thus LIKWID reads correctly the right hardware counters.
