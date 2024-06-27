# LeMonADE_StiffCyclic
Extension of LeMonADE package for the simulaton of stiff, cyclic polymers.

## projects/createTwoConcatenatedRings/
Creation of two concatenated ring polymers of arbitary length (N,M >= 8) in a non-periodic simulation box.  

* To create two concatenated rings with M=50 monomers and N=100 monomers in V=256³ simulation box:
```
    ./CreateTwoConcatenatedRings -o TwoConcatenatedRings_M50_N100_NoPerXYZ256.bfm -m 50 -n 100 -b 256
```
