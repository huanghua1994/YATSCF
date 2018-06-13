#ifndef _YATSCF_BUILD_DENSITY_H_
#define _YATSCF_BUILD_DENSITY_H_

#include "TinySCF.h"

// Build density matrix using eigen decomposition (diagonalization)
void TinySCF_build_DenMat(TinySCF_t TinySCF);

// Build density matrix using Canonical Purification
#define PURIF_TOL      1e-11
#define MAX_PURIF_ITER 200
void TinySCF_build_DenMat_Purif(TinySCF_t TinySCF, int *purif_iter);

// Build density matrix using Second-order spectral projection (SP2)
#define SP2_TOL        1e-11
#define MAX_SP2_ITER   200
void TinySCF_build_DenMat_SP2(TinySCF_t TinySCF, int *SP2_iter);

#endif
