#ifndef _YATSCF_BUILD_DENSITY_H_
#define _YATSCF_BUILD_DENSITY_H_

#include "TinySCF.h"

// Build density matrix using eigen decomposition (diagonalization)
void TinySCF_build_DenMat(TinySCF_t TinySCF);

// Build density matrix using Canonical Purification
#define PURIF_TOL      1e-11
#define MAX_PURIF_ITER 200
void TinySCF_build_DenMat_Canonical(TinySCF_t TinySCF, int *purif_iter);

// Build density matrix using Second-order spectral projection (SP2)
#define SP2_TOL        1e-11
#define MAX_SP2_ITER   200
void TinySCF_build_DenMat_SP2(TinySCF_t TinySCF, int *SP2_iter);

// Build density matrix using Stable, Scaled Newton-Schulz method (SSNS)
#define SSNS_TOL       1e-11
#define MAX_SSNS_ITER  200
void TinySCF_build_DenMat_SSNS(TinySCF_t TinySCF, int *SSNS_iter);

// Build density matrix using McWeeny Purification
#define MCWEENY_TOL       1e-11
#define MAX_MCWEENY_ITER  200
void TinySCF_build_DenMat_McWeeny(TinySCF_t TinySCF, int *McWeeny_iter);

#endif
