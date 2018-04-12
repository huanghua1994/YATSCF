#ifndef _YATSCF_ACCUM_FOCK_H_
#define _YATSCF_ACCUM_FOCK_H_

#include "TinySCF.h"

void Accum_Fock(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

#endif