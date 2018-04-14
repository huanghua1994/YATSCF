/*
 * Copyright (c) 2013-2018 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#ifndef __CINT_SIMINT_H__
#define __CINT_SIMINT_H__

struct SIMINT
{
    int nthreads;
    int max_am;
    int workmem_per_thread;
    int outmem_per_thread;
    double *workbuf;
    double *outbuf;

    int nshells;
    struct simint_shell *shells;
    struct simint_multi_shellpair *shellpairs;

    int    screen_method;
    double screen_tol;

    // For timer, only master thread will write to these
    double ostei_actual, ostei_setup, fock_update_F;

    // For statistic
    double *num_multi_shellpairs, *sum_nprim;
    double *num_screened_prim, *num_unscreened_prim, *num_screened_vec, *num_unscreened_vec;
};

typedef struct SIMINT *SIMINT_t;

#ifdef __cplusplus
extern "C" {
#endif

CIntStatus_t CInt_createSIMINT(BasisSet_t basis, SIMINT_t *simint, int nthreads);

CIntStatus_t CInt_destroySIMINT(SIMINT_t simint, int show_stat);

CIntStatus_t
CInt_computeShellQuartet_SIMINT(SIMINT_t simint, int tid,
                                int A, int B, int C, int D,
                                double **integrals, int *nints);

CIntStatus_t
CInt_computePairOvl_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                           int A, int B,
                           double **integrals, int *nints);

CIntStatus_t
CInt_computePairCoreH_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                           int A, int B,
                           double **integrals, int *nints);


// The following 2 constants are corresponding to SIMINT_OSTEI_MAXAM
// and SIMINT_NSHELL_SIMD in Simint. I cannot include <simint/simint.h>
// here, so I just update the values manually. This problem should be 
// solved later.
#define _SIMINT_OSTEI_MAXAM 4
#define _SIMINT_NSHELL_SIMD 16

#define _SIMINT_AM_PAIRS (((_SIMINT_OSTEI_MAXAM) + 1) * ((_SIMINT_OSTEI_MAXAM) + 1))

void CInt_SIMINT_addupdateFtimer(SIMINT_t simint, double sec);

int  CInt_SIMINT_getShellpairAMIndex(SIMINT_t simint, int P, int Q);

void CInt_SIMINT_createThreadMultishellpair(void **thread_multi_shellpair);

void CInt_SIMINT_freeThreadMultishellpair(void **thread_multi_shellpair);

CIntStatus_t 
CInt_computeShellQuartetBatch_SIMINT(
    SIMINT_t simint, int tid,
    int M, int N, int *P_list, int *Q_list,
    int npair, double **thread_batch_integrals, int *thread_batch_nints,
    void **thread_multi_shellpairs
);

#ifdef __cplusplus
}
#endif

#endif