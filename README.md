# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## Installation

* Clone and Install LeMonADE v2.2.2, see [LeMonADE v2.2.2](https://github.com/LeMonADE-project/LeMonADE/tree/v2.2.2)
* Install cmake (minimum version 3.1)
* Just do for standard compilation of the ELMA-OscillatoryForce project:

```sh
    # generates the projects
    mkdir build
    cd build
    cmake -DLEMONADE_DIR=/path/to/LeMonADE-library/ ..
    make
```

## projects/AnalyzerScatteringSingleObject/

Calculates the scattering intensity S(q) with respect to the absolute magnitute of the scattering wave vector q.  

```
    ./AnalyzerScatteringSingleObject [-f filename(=test.bfm)] [-e evaluation_time(=0)] [-s skipframe(=0)]
```

## projects/createTwoConcatenatedRings/
Creation of two concatenated ring polymers of arbitary length (N,M >= 8) in a non-periodic simulation box.  

* To create two concatenated rings with M=50 monomers and N=100 monomers in V=256³ simulation box:
```
    ./CreateTwoConcatenatedRings -o TwoConcatenatedRings_M50_N100_NoPerXYZ256.bfm -m 50 -n 100 -b 256
```
