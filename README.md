# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Instructions for Generating Simulation Executables
1. Install cmake (minimum version 3.1).
2. Clone and Install LeMonADE v2.2.2, see [LeMonADE v2.2.2](https://github.com/LeMonADE-project/LeMonADE/tree/v2.2.2).
3. Compile LeMonADE_StiffCyclic programs.
```
./results/00_programs/bash/build_LeMonADE.sh -p <path/to/lemonade>
```
4. Generate initial polymer configuration and settings.
```
./build/bin/generatePolymerBFM -h

Usage: ./generatePolymerBFM << OPTIONS >> 
[-o filenameOutput] 
[-n number of monomers in ring / chain] 
[-m number of rings / chains] 
[-r generate ring (otherwise generate chain)] 
[-k bending potential strength (otherwise no bending potential)] 
[-f constant force that molecules experience (otherwise no force is applied)] 
[-v string with three integers xyz denoting the force orientation in each dimension (default is 111)] 
[-b box size]
[-c use cosine angle potential for bending potential (default is cosine square angle potential)]
[-p force oscillation period (default is 1)]
[-a force oscillation amplitude (1)]

```
5. Simulate Polymer system.
```
./build/bin/simulatePolymerBFM -h

Usage: ./simulatePolymerBFM << OPTIONS >> 
[-i load file] 
[-o output file] 
[-n number of total Monte Carlo steps] 
[-s save frequency (in Monte Carlo steps)]
[-e equilibriation time]
[-a add end-to-end analyzer]
[-b add hysteresis analyzer]
[-c add bond-bond correlation analyzer]
[-d add bond vector distribution analyzer]
[-g add radius of gyration analyzer]

```
## Submitting batch jobs
>[!NOTE]
>These programs were designed to execute on the IPFDD remote computing cluster system.
1. Navtigate to `./results`
2. Generate instructions and directory hirearchy for parameterzed simulations in `./01_raw_data`.
	a.  *bendingPARM*: Generate simulations in which the bending potential strength is parameterized. 
	b.  *forceExtension*: Generate force extension simulations in which stiff polymers under a constant force. 
	c.  *hysteresis*: Generate hysteresis simulations in which polymers of various stiffness and topologies are exposed to an oscillating sinusodal force.
3. Run batch simulations.
	*  Run locally.
	*  Submit to computing cluster.
4. Analyze simulation results.

