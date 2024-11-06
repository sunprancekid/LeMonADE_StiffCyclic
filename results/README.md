# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Instructions for Generating, Submitting, and Analyzing Batch Simulations
>[!NOTE]
>These programs were designed to execute on the IPFDD remote computing cluster system.

Generally speaking, simulations are generated in batch according to the proceedure outlined below.
1. *Generate Simulation Parameters and Directories*: There are three unqiue modules which generate simulations according to a unique set of parameters that are varied (outline below). Each module has a unique job code (`JOBID`). By calling, `./00_programs/bash/gen_${JOBID}.sh -g` simulations and their parameters are generated and stored within `./01_raw_data/${JOBID}/`. Each simulation has a unique set of parameters, a unique simulation code (`SIMID`), and a unique directory (`SIMDIR`). All of the parameters associated with one module / job are stored within `./01_raw_data/${JOBID}/${JOBID}.csv`.
2. *Run Simulations*: Run simulations in batch, either locally on you machine (not recommended) or on a remote computing cluster.
3. *Compile and Analyze Simulation Results*: Analyze simulation results.

There are three different modules which can be used to generate sets of simulations in which certain sets of parameters are varied in order to explore their effect of simulated polymer systems.
###bendingPARM
Generate simulations in which the bending potential strength is parameterized. Explore the effect of the polymer stiffness of ideal and real polymer chains.

###forceExtension
Generate force extension simulations in which stiff polymers under a constant force. 

###hysteresis
Generate hysteresis simulations in which polymers of various stiffness and topologies are exposed to an oscillating sinusodal force.
