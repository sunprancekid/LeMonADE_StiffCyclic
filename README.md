# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Instructions for Generating Simulation Executables
1. Install CMake.
2. Install [LeMonADE Package](https://github.com/LeMonADE-project/LeMonADE).
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
