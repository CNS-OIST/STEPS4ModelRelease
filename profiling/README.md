Here, we describe how to perform the strong scaling & Roofline analysis found in STEPS 4 Release paper. To simply re-generate the graphs of the paper just run the `STEPS_PerfGraphs.ipynb` jupyter notebook (all data is baked in).

## How to perform a Strong Scaling analysis

There are 3 main case studies, namely SimpleModel, caBurstBackground, and caBurstFull, for which we are interested in performing a strong scaling analysis. To perform a strong scaling analysis, the problem size remains fixed while the number of cores is steadily increased. The steps to do so are presented per case study below:

1. SimpleModel (see corresponding folder `../SimpleModel/profiling`): Modify accordingly and run the `CI_Perf.ipynb` jupyter notebook, which does the setup of the environment, the job submission, and the graph generation (adds also a comparison with a reference STEPS version, in our case it can be the version of STEPS 4 Release paper -found in `perfRes_STEPS4.0Paper` folder-). This notebook is integrated in our Continuous Integration (CI) pipeline and is supposed to be executed by a GitLab runner, i.e., from an already allocated job. Therefore, some modifications must be done to adjust it to your environment.

2. caBurstBackground (see corresponding folder `../caBurst/background/profiling`): Exact same steps as in **SimpleModel**.

3. caBurstFull (see corresponding folder `../caBurst/full/profiling`): Exact same steps as in **SimpleModel**

### Strong Scaling analysis : Caliper Instrumentation

The same notebook (`CI_Perf.ipynb`) performs the Caliper Instrumentation as well. Keep in mind that STEPS must be compiled with Caliper support `ON` (see CMake variable `STEPS_USE_CALIPER_PROFILING`).

## How to perform a Roofline analysis

We have chosen to perform the Roofline analysis for the caBurstFull case study given that it combines all the main computational kernels, i.e., SSA operator, Diffusion operator, and Efield operator. In the `../caBurst/full/profiling/Roofline_Analysis/` folder, one can find the job submission script (always single node analysis):

* Submit multiple jobs in the queue of a supercomputer (in our case is BB5 of the Blue Brain Project) with the help of `likwid.batch` script. The python script to run the caBurstFull model is the same as the one of strong scaling (can be found in `../caBurst/full/profiling/`):
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
1. STEPS must be compiled with LIKWID support `ON` (see CMake variable `STEPS_USE_LIKWID_PROFILING`).

2. It is very crucial to use `export STEPS_INSTRUMENTOR_MPI_BARRIER="before;after"` to make sure that we profile with LIKWID the annotated regions correctly. This happens in `likwid.batch` job submission script. The MPI barrier guarantees that all the processes traverse at the same time the annotated computational kernels, and thus LIKWID reads correctly the right hardware counters.
