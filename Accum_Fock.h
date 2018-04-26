#ifndef _YATSCF_ACCUM_FOCK_H_
#define _YATSCF_ACCUM_FOCK_H_

#include "TinySCF.h"

#define ACCUM_FOCK_IN_PARAM TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, \
							double *ERI, int load_MN, int load_P, int write_MN, int write_P

void Accum_Fock(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_dimQ1(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_dimQ3(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_dimQ6(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_dimQ10(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_dimQ15(ACCUM_FOCK_IN_PARAM);

void Accum_Fock_1111(ACCUM_FOCK_IN_PARAM);

#endif