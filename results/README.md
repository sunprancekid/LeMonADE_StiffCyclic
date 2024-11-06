# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Instructions for Generating, Submitting, and Analyzing Batch Simulations

Generally speaking, simulations are generated in batch from the command line terminal (CLT) according to the proceedure outlined below.
1. *Generate Simulation Parameters and Directories*: There are three unqiue modules which generate simulations according to a unique set of parameters that are varied (outline below). Each module has a unique job code (`JOBID`). By calling, `./00_programs/bash/gen_${JOBID}.sh -g` simulations and their parameters are generated and stored within `./01_raw_data/${JOBID}/`. Each simulation has a unique set of parameters, a unique simulation code (`SIMID`), and a unique directory (`SIMDIR`). All of the parameters associated with one module / job are stored within `./01_raw_data/${JOBID}/${JOBID}.csv`.
2. *Run Simulations*: Run simulations in batch, either locally on you machine by calling `./00_programs/bash/submit_batch_jobs.sh -j ${JOBID}` (not recommended) or on a remote computing cluster `./00_programs/bash/submit_batch_jobs.sh -s -j ${JOBID}`. The flag `-l` can be used to specify a remote cluster from the default.
>[!NOTE]
>These scripts were designed to submit simulations on the IPFDD remote computing cluster systems.
>Furthermore, the `LeMonADE_StiffCyclic` should also be installed on comupting cluster in question.
```
USAGE: ./submit_batch_job.sh << FLAGS >>

 ## SCRIPT PROTOCOL ##
 -h           | display script options, exit 0.
 -v           | execute script verbosely.
 -s           | submit job to SLURM via remote linux cluster (otherwise run locally in serial).
 -t           | submit one job as test in order to make sure everthing works properly.

 ## SCRIPT PARAMETERS ##
 -j << ARG >> | MANDATORY: specify job title.
 -p << ARG >> | specify path to simulation job directory hirearchy (default is 01_raw_data/${JOB}).
 -f << ARG >> | file which contains simulation parameters (default is 01_raw_data/${JOB}/${JOB}.csv)
 -l << ARG >> | specify remote linux cluster (default is gandalf).
 -e << ARG >> | specify path to LeMonADE executables (default is 00_programs/build/bin/).
 -c << ARG >> | check for file; if this file exists in the simulation directory, do not execute.

```

3. *Compile and Analyze Simulation Results*: Analyze simulation results.

There are three different modules which can be used to generate sets of simulations in which certain sets of parameters are varied in order to explore their effect of simulated polymer systems.

### bendingPARM
Generate simulations in which the bending potential strength is parameterized. Explore the effect of the polymer stiffness of ideal and real polymer chains.

### forceExtension
Generate force extension simulations in which stiff polymers under a constant force. 

### `JOBID=hysteresis`
Generate hysteresis simulations in which polymers of various stiffness and topologies are exposed to an oscillating sinusodal force.
