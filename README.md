# YATSCF: Yet Another Tiny SCF

My final project for CSE 8803 Spring 2018, a HF-SCF program:

* Written in C, Parallelized with OpenMP in single node ("Tiny");
* Used [Simint](https://github.com/simint-chem/simint-generator) for ERIs;
* Accelerate SCF convergence with SAD initial guess and DIIS.

In short, this program can be viewed as a tiny, more "formula translation" version of [GTFock](https://github.com/gtfock-chem/gtfock).

Current status:

![](https://img.shields.io/badge/Correct%20results-yes-brightgreen.svg)
![](https://img.shields.io/badge/SAD%20init%20guess-not%20ready-red.svg)
![](https://img.shields.io/badge/DIIS-slow%20conv.-yellow.svg)
![](https://img.shields.io/badge/ERI%20batching-not%20ready-red.svg)
![](https://img.shields.io/badge/OpenMP%20parallelization-not%20ready-red.svg)
![](https://img.shields.io/badge/Fock%20accum.%20opt.-not%20ready-red.svg)