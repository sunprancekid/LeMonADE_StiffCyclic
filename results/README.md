# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Instructions for Generating, Submissing, and Analyzing Batch Simulations
>[!NOTE]
>These programs were designed to execute on the IPFDD remote computing cluster system.
Generally speaking, simulations can be generated in batch according to the proceedure outlined below.
1. Generate instructions and directory hirearchy for parameterzed simulations in `./01_raw_data`. There are three
	*  *bendingPARM*: Generate simulations in which the bending potential strength is parameterized. 
	*  *forceExtension*: Generate force extension simulations in which stiff polymers under a constant force. 
	*  *hysteresis*: Generate hysteresis simulations in which polymers of various stiffness and topologies are exposed to an oscillating sinusodal force.
2. Run simulations in batch, either locally on you machine (not recommended) or on a remote computing cluster.
3. Analyze simulation results.

There are three different modules which can be used to generate sets of simulations in which certain sets of parameters are varied in order to explore their effect of simulated polymer systems.
