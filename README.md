# YATSCF: Yet Another Tiny SCF

My final project for CSE 8803 Spring 2018, a HF-SCF program:

* Written in C99, parallelized with OpenMP in single node ("Tiny");
* Used [Simint](https://github.com/simint-chem/simint-generator) for ERIs;
* Used SAD initial guess and DIIS(Pulay mixing) to accelerate SCF convergence.

In short, this program can be viewed as a tiny, more "formula translation" version of [GTFock](https://github.com/gtfock-chem/gtfock).

Current status:

![](https://img.shields.io/badge/Same%20conv.%20energy%20as%20GTFock-yes-brightgreen.svg)

![](https://img.shields.io/badge/SAD%20initial%20guess-ready-brightgreen.svg)

![](https://img.shields.io/badge/DIIS%20acceleration-ready-brightgreen.svg)

![](https://img.shields.io/badge/OpenMP%20parallelization-ready-brightgreen.svg)

![](https://img.shields.io/badge/ERI%20batching-ready-brightgreen.svg)

![](https://img.shields.io/badge/Fock%20accum.%20opt.-ready-brightgreen.svg)

![](https://img.shields.io/badge/Repack%20libcint-no-red.svg)
