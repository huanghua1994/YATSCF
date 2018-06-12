#ifndef _YATSCF_BUILD_DENSITY_H_
#define _YATSCF_BUILD_DENSITY_H_

#include "TinySCF.h"

// Build density matrix using eigen decomposition (diagonalization)
void TinySCF_build_DenMat(TinySCF_t TinySCF);

// Build density matrix using Canonical Purification
#define MAX_PURIF_ITER 200
void TinySCF_build_DenMat_Purif(TinySCF_t TinySCF, int *purif_iter);

#endif
