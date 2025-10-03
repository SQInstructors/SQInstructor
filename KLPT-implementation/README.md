# KLPT-based implementation


This code is mostly meant for playing around with the KLPT variants our paper proposes for Scalar and Borel level structures. 


## Origin of code

Most parts of the code, in particular the ideal-to-isogeny related parts and many KLPT subroutines were copied from [LearningToSQI](https://github.com/LearningToSQI/SQISign-SageMath) by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer and Giacomo Pope.
This includes many subfunctions of KLPT which were partially copied to new files and edited.
In summary, all code in the mytest_...sage files, and parts of the content of ScalarKLPT.py, BorelKLPT.py, are from us, most of the rest from LearningToSQI.

## Requirements

- A recent version of Sagemath (10.5 is fine)
- Fpylll  (0.6.1 s fine)

## Structure

Short overview over the many files here.

### Infrastructure

Many smaller and larger helper files, mostly taken from LearningToSQI
- setup, parameters: Signature and KLPT parameter settings which are used in most files
- helpers, lattices, ideals, utilities, rerandomization: Small functions, rather on the quaternion side
- mitm, isogenies, deuring, torsion: ideal-to-isogeny tools
- pari_interface: Low-level plumbing


### KLPT variants

- ScalarKLPT.py Implements the (generalized) algorithm for KLPT respecting scalar level structures
- BorelKLPT.py Implements the most efficient version of KLPT respecting a borel level structure
None of these implementation has well-distributed outputs (or even just fixed output norm), and both fail from time to time, for example on inputs with large forth minimum. 
To summarize, both of them only support power-of-two output norms, odd power-of-prime values of the torsion integer Ntor, and always use the special extremal order O0 (of basis1,i,(j+i)/2,(1+ij)/2)for solving. Also inputs must have norm coprime to Ntor.


### 1d SQIsign implementation

- sqitorsion.py: Implementes the efficient a signature scheme around KLPT. So far only reliably works for Borel level structures.
### Tests

- mytest_klpt.sage: Runs a simple, quaternion-only test of both KLPT variants
- mytest_isogenies.sage: Runs a simple test of torsion mapping, on ScalarKLPT and BorelKLPT
- mytest_sqitorsion.sage: Runs a sqitorsion key generation, signing and verification, using Borel level structures
