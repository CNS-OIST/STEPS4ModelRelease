# STEPS4ModelRelease
### Introduction

In this repository we present the models described in the STEPS 4 paper. In order to run the models you need clone and build STEPS 4 from: https://github.com/CNS-OIST/STEPS. 

##### Installation of STEPS

Without loss of generality, let us assume that you installed STEPS in `~/STEPS` and the build folder in  `~/STEPS/build`. In order to make steps discoverable by python we need to add `~/STEPS/build/lib` to the pythonpath using the following bash command:

```bash
export PYTHONPATH=~/STEPS/build/lib:$PYTHONPATH
```

After, you are ready to run the models. 

##### Running the models

For more information on how to run a specific model we suggest to check the specific README.md in the model folder. In every case, at the end you get a file `res*.txt` in the respective `results` folder. It contains the raw traces that can be statistically analyzed. 

##### Validations

To analyse the raw traces we suggest to download https://github.com/CNS-OIST/STEPS_Validation and head toward the folder `postproc`. In there you can copy the traces and perform some statistical analyses. Check the relative README and the examples for more information. 

##### Profiling/ Performance

The measurement of STEPS performance is achieved through suitable instrumentation. To simply generate all the graphs found in STEPS 4 Release paper, locate the `STEPS_PerfGraphs.ipynb` jupyter notebook in the `./profiling` folder, and run all the cells of the notebook. All data needed is already baked in the notebook. However, the interested reader/scientist could follow the instructions found in `./profiling/README.md` and perform the strong scaling & Roofline analysis to re-generate similar data, using for example a different supercomputing facility.
