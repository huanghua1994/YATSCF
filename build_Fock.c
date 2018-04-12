#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "utils.h"
#include "TinySCF.h"
#include "build_Fock.h"
#include "shell_quartet_list.h"

static inline void atomic_add_f64(volatile double *global_value, double addend)
{
	uint64_t expected_value, new_value;
	do {
		double old_value = *global_value;
		expected_value = _castf64_u64(old_value);
		new_value = _castf64_u64(old_value + addend);
	} while (!__sync_bool_compare_and_swap((volatile uint64_t*)global_value, expected_value, new_value));
}

static inline void unique_integral_coef(int M, int N, int P, int Q, double *coef)
{
	int flag1 = (M == N) ? 0 : 1;
	int flag2 = (P == Q) ? 0 : 1;
	int flag3 = ((M == P) && (N == Q)) ? 0 : 1;
	int flag4 = ((flag1 == 1) && (flag2 == 1)) ? 1 : 0;
	int flag5 = ((flag1 == 1) && (flag3 == 1)) ? 1 : 0;
	int flag6 = ((flag2 == 1) && (flag3 == 1)) ? 1 : 0;
	int flag7 = ((flag4 == 1) && (flag3 == 1)) ? 1 : 0;
	coef[0] = 2.0 * (1.0   + flag1 + flag2 + flag4);  // for J_MN
	coef[1] = 2.0 * (flag3 + flag5 + flag6 + flag7);  // for J_PQ
	coef[2] = 1.0   + flag3;  // for K_MP
	coef[3] = flag1 + flag5;  // for K_NP
	coef[4] = flag2 + flag6;  // for K_MQ
	coef[5] = flag4 + flag7;  // for K_NQ
}

void Accum_Fock(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI
)
{
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = TinySCF->shell_bf_num[Q];
	int idxM = TinySCF->shell_bf_sind[M];
	int idxN = TinySCF->shell_bf_sind[N];
	int idxP = TinySCF->shell_bf_sind[P];
	int idxQ = TinySCF->shell_bf_sind[Q];
	
	double *J_MN = TinySCF->J_mat + idxM * nbf + idxN;
	double *J_PQ = TinySCF->J_mat + idxP * nbf + idxQ;
	double *K_MP = TinySCF->K_mat + idxM * nbf + idxP;
	double *K_NP = TinySCF->K_mat + idxN * nbf + idxP;
	double *K_MQ = TinySCF->K_mat + idxM * nbf + idxQ;
	double *K_NQ = TinySCF->K_mat + idxN * nbf + idxQ;
	
	double *D_MN = TinySCF->D_mat + idxM * nbf + idxN;
	double *D_PQ = TinySCF->D_mat + idxP * nbf + idxQ;
	double *D_MP = TinySCF->D_mat + idxM * nbf + idxP;
	double *D_NP = TinySCF->D_mat + idxN * nbf + idxP;
	double *D_MQ = TinySCF->D_mat + idxM * nbf + idxQ;
	double *D_NQ = TinySCF->D_mat + idxN * nbf + idxQ;
	
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	// TODO: leading dimension of imn, inp, imp, ipq_base, imq_base, inq_base
	// will no longer be nbf when using thread-private buffer
	for (int iM = 0; iM < dimM; iM++) 
	{
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iMN = iM * nbf + iN;
			int iMQ_base = iM * nbf;
			int iNQ_base = iN * nbf;
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iMP = iM * nbf + iP;
				int iNP = iN * nbf + iP;
				int iPQ_base = iP * nbf;
				int Ibase = dimQ * (iP + dimP * (iN + dimN * iM));
				double k_MP = 0.0, k_NP = 0.0;
				// dimQ is small, vectorizing short loops may hurt performance since
				// it needs horizon reduction after the loop
				for (int iQ = 0; iQ < dimQ; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += coef[0] * D_PQ[iPQ_base + iQ] * I;
					k_MP -= coef[2] * D_NQ[iNQ_base + iQ] * I;
					k_NP -= coef[3] * D_MQ[iMQ_base + iQ] * I;
					
					atomic_add_f64(&J_PQ[iPQ_base + iQ], coef[1] * D_MN[iMN] * I);
					atomic_add_f64(&K_MQ[iMQ_base + iQ], -coef[4] * D_NP[iNP] * I);
					atomic_add_f64(&K_NQ[iNQ_base + iQ], -coef[5] * D_MP[iMP] * I);
				}
				atomic_add_f64(&K_MP[iMP], k_MP);
				atomic_add_f64(&K_NP[iNP], k_NP);
			} // for (int iM = 0; iM < dimM; iM++) 
			atomic_add_f64(&J_MN[iMN], j_MN);
		} // for (int iQ = 0; iQ < dimQ; iQ++) 
	} // for (int iN = 0; iN < dimN; iN++)
}

void Accum_Fock_with_KetshellpairList(
	TinySCF_t TinySCF, int tid, int M, int N, 
	int *P_list, int *Q_list, int npairs, 
	double *thread_eris, int thread_nints
)
{
	for (int ipair = 0; ipair < npairs; ipair++)
	{
		Accum_Fock(
			TinySCF, tid, M, N, P_list[ipair], Q_list[ipair], 
			thread_eris + ipair * thread_nints
		);
	}
}

// F = H_core + (J + J^T) / 2 + (K + K^T) / 2;
void TinySCF_HJKmat_to_Fmat(double *Hcore_mat, double *J_mat, double *K_mat, double *F_mat, int nbf)
{	
	for (int irow = 0; irow < nbf; irow++)
	{
		int idx = irow * nbf + irow;
		F_mat[idx] = Hcore_mat[idx] + J_mat[idx] + K_mat[idx];
		for (int icol = irow + 1; icol < nbf; icol++)
		{
			int idx1 = irow * nbf + icol;
			int idx2 = icol * nbf + irow;
			double Jval = (J_mat[idx1] + J_mat[idx2]) * 0.5;
			double Kval = (K_mat[idx1] + K_mat[idx2]) * 0.5;
			F_mat[idx1] = Hcore_mat[idx1] + Jval + Kval;
			F_mat[idx2] = Hcore_mat[idx2] + Jval + Kval;
		}
	}
}

void TinySCF_build_FockMat(TinySCF_t TinySCF)
{
	// Copy some parameters out, I don't want to see so many "TinySCF->"
	int nshells        = TinySCF->nshells;
	int num_uniq_sp    = TinySCF->num_uniq_sp;
	double scrtol2     = TinySCF->shell_scrtol2;
	double *sp_scrval  = TinySCF->sp_scrval;
	int *shell_bf_num  = TinySCF->shell_bf_num;
	int *shell_bf_sind = TinySCF->shell_bf_sind;
	int *uniq_sp_lid   = TinySCF->uniq_sp_lid;
	int *uniq_sp_rid   = TinySCF->uniq_sp_rid;
	int num_bas_func   = TinySCF->nbasfuncs;
	SIMINT_t simint    = TinySCF->simint;
	double *J_mat      = TinySCF->J_mat;
	double *K_mat      = TinySCF->K_mat;
	double *F_mat      = TinySCF->F_mat;
	double *D_mat      = TinySCF->D_mat;
	double *Hcore_mat  = TinySCF->Hcore_mat;
	
	memset(J_mat, 0, DBL_SIZE * TinySCF->mat_size);
	memset(K_mat, 0, DBL_SIZE * TinySCF->mat_size);
	
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		
		// Create ERI batching auxiliary data structures
		// Ket-side shell pair lists that needs to be computed
		ThreadKetShellpairLists_t thread_ksp_lists;
		create_ThreadKetShellpairLists(&thread_ksp_lists);
		// Simint multi_shellpair buffer for batched ERI computation
        void *thread_multi_shellpair;
        CInt_SIMINT_createThreadMultishellpair(&thread_multi_shellpair);
		
		#pragma omp for schedule(dynamic)
		for (int MN = 0; MN < num_uniq_sp; MN++)
		{
			int M = uniq_sp_lid[MN];
			int N = uniq_sp_rid[MN];
			double scrval1 = sp_scrval[M * nshells + N];
			for (int PQ = 0; PQ < num_uniq_sp; PQ++)
			{
				int P = uniq_sp_lid[PQ];
				int Q = uniq_sp_rid[PQ];
				double scrval2 = sp_scrval[P * nshells + Q];
				
				// Symmetric uniqueness check, from GTFock
				if ((M > P && (M + P) % 2 == 1) || 
					(M < P && (M + P) % 2 == 0))
				continue; 
				
				if ((M == P) &&	((N > Q && (N + Q) % 2 == 1) ||
					(N < Q && (N + Q) % 2 == 0)))
				continue;
				
				// Shell screening 
				if (fabs(scrval1 * scrval2) <= scrtol2) continue;
				
				// Push ket-side shell pair to corresponding list
				int ket_id = CInt_SIMINT_getShellpairAMIndex(simint, P, Q);
				KetShellpairList_t target_shellpair_list = &thread_ksp_lists->ket_shellpair_lists[ket_id];
				add_shellpair_to_KetShellPairList(target_shellpair_list, P, Q);
				
				// If the ket-side shell pair list we just used is full, handle it
				if (target_shellpair_list->npairs == MAX_LIST_SIZE)
				{
					double *thread_batch_eris;
					int thread_nints;
					
					// Compute batched ERIs
					CInt_computeShellQuartetBatch_SIMINT(
						simint, tid, M, N,
						target_shellpair_list->P_list,
						target_shellpair_list->Q_list,
						target_shellpair_list->npairs,
						&thread_batch_eris, &thread_nints, 
						&thread_multi_shellpair
					);
					
					// Accumulate ERI results to global matrices
					if (thread_nints > 0)
					{
						double st = get_wtime_sec();
						Accum_Fock_with_KetshellpairList(
							TinySCF, tid, M, N, 
							target_shellpair_list->P_list,
							target_shellpair_list->Q_list,
							target_shellpair_list->npairs,
							thread_batch_eris, thread_nints
						);
						double et = get_wtime_sec();
						if (tid == 0) CInt_SIMINT_addupdateFtimer(simint, et - st);
					}
					
					// Reset the computed ket-side shell pair list
					target_shellpair_list->npairs = 0;
				}
			}
			
			// Handles all non-empty ket-side shell pair lists
			for (int ket_id = 0; ket_id < MAX_AM_PAIRS; ket_id++)
			{
				KetShellpairList_t target_shellpair_list = &thread_ksp_lists->ket_shellpair_lists[ket_id];
				
				if (target_shellpair_list->npairs > 0)
				{
					double *thread_batch_eris;
					int thread_nints;
					
					// Compute batched ERIs
					CInt_computeShellQuartetBatch_SIMINT(
						simint, tid, M, N,
						target_shellpair_list->P_list,
						target_shellpair_list->Q_list,
						target_shellpair_list->npairs,
						&thread_batch_eris, &thread_nints, 
						&thread_multi_shellpair
					);
					
					// Accumulate ERI results to global matrices
					if (thread_nints > 0)
					{
						double st = get_wtime_sec();
						Accum_Fock_with_KetshellpairList(
							TinySCF, tid, M, N, 
							target_shellpair_list->P_list,
							target_shellpair_list->Q_list,
							target_shellpair_list->npairs,
							thread_batch_eris, thread_nints
						);
						double et = get_wtime_sec();
						if (tid == 0) CInt_SIMINT_addupdateFtimer(simint, et - st);
					}
					
					// Reset the computed ket-side shell pair list
					target_shellpair_list->npairs = 0;
				}
			}
		}
		
		// Free batch ERI auxiliary data structures
		CInt_SIMINT_freeThreadMultishellpair(&thread_multi_shellpair);
		free_ThreadKetShellpairLists(thread_ksp_lists);
	}
	
	TinySCF_HJKmat_to_Fmat(Hcore_mat, J_mat, K_mat, F_mat, num_bas_func);
}
