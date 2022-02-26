# STEPS4 simulation models for submission to Frontiers in Neuroinformatics ("STEPS 4.0: Fast and memory-efficient molecular simulations of neurons at the nanoscale")
### Introduction

In this repository we present the models described in the manuscript submitted to Frontiers in Neuroinformatics titled "STEPS 4.0: Fast and memory-efficient molecular simulations of neurons at the nanoscale". 

Please note that the solution described in the manuscript focusing on large scale cpu-based high performance computing clusters, thus some of the models may not be suitable for regular desktop machines. We also provide the data we gathered in our study for readers who are interested in statistical analysis of the existing results.   

In order to run the models an installation of STEPS 4.0 is required. Please clone and build STEPS 4.0 from: https://github.com/CNS-OIST/STEPS. 

##### Installation of STEPS

Without loss of generality, let us assume that you installed STEPS in `/path/to/STEPS`, see [STEPS installation documentation](https://github.com/CNS-OIST/STEPS/#installation-from-source-code) for more details. In order to make the simulator discoverable by Python, it is required to provide its location using the following Bash command:

```bash
export PYTHONPATH="/path/to/STEPS:$PYTHONPATH"
```

After, you are ready to run the models. Additionally, you can:

```bash
export STEPS_INSTALL_DIR=/path/to/STEPS
```

and, use it to locate the installed STEPS.

##### Running the models

For more information on how to run a specific model we suggest to check the specific README.md in the model folder. In every case, at the end a file `res*.txt` is generated in the respective `results` folder. It contains the raw traces that can be statistically analyzed. 

##### Validations

To analyse the raw traces we suggest to download https://github.com/CNS-OIST/STEPS_Validation and head toward the folder `postproc`. In there you can copy the traces and perform some statistical analyses. Check the relative README and the examples for more information. 

##### Profiling/ Performance

The measurement of STEPS performance is achieved through suitable instrumentation. To simply generate all the graphs found in the submitted manuscript, locate the `STEPS_PerfGraphs.ipynb` jupyter notebook in the `./profiling` folder, and run all the cells of the notebook. All data needed is already baked in the notebook. However, the interested reader/scientist could follow the instructions found in `./profiling/README.md` and perform the strong scaling & Roofline analysis to re-generate similar data, using for example a different supercomputing facility.
