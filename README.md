# YATSCF: Yet Another Tiny SCF

My final project for CSE 8803 Spring 2018, a HF-SCF program:

* Written in C99, parallelized with OpenMP in single node ("Tiny");
* Used [Simint](https://github.com/simint-chem/simint-generator) for ERIs;
* Accelerate SCF convergence with SAD initial guess and DIIS.

In short, this program can be viewed as a tiny, more "formula translation" version of [GTFock](https://github.com/gtfock-chem/gtfock).

Current status:

![](https://img.shields.io/badge/Same%20conv.%20energy%20as%20GTFock-yes-brightgreen.svg)

![](https://img.shields.io/badge/SAD%20initial%20guess-ready-brightgreen.svg)

![](https://img.shields.io/badge/DIIS%20acceleration-problematic-yellow.svg)

![](https://img.shields.io/badge/OpenMP%20parallelization-ready-brightgreen.svg)

![](https://img.shields.io/badge/ERI%20batching-ready-brightgreen.svg)

![](https://img.shields.io/badge/Fock%20accum.%20opt.-partly%20ready-green.svg)
