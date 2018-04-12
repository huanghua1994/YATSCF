#ifndef _YATSCF_ACCUM_FOCK_H_
#define _YATSCF_ACCUM_FOCK_H_

#include "TinySCF.h"

void Accum_Fock(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_dimQ1(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_dimQ3(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_dimQ6(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_dimQ10(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_dimQ15(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

void Accum_Fock_1111(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
);

#endif