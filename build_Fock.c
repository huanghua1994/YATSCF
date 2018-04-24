#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "utils.h"
#include "TinySCF.h"
#include "build_Fock.h"
#include "shell_quartet_list.h"
#include "Accum_Fock.h"
#include "CMS.h"

#define ACCUM_FOCK_PARAM 	TinySCF, tid, M, N, P_list[ipair], Q_list[ipair], \
							thread_eris + ipair * thread_nints, \
							load_MN, load_P, write_MN, write_P

void Accum_Fock_with_KetshellpairList(
	TinySCF_t TinySCF, int tid, int M, int N, 
	int *P_list, int *Q_list, int npairs, 
	double *thread_eris, int thread_nints
)
{
	int load_MN, load_P, write_MN, write_P, prev_P = -1;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	for (int ipair = 0; ipair < npairs; ipair++)
	{
		load_MN = (ipair == 0) ? 1 : 0;
		load_P  = (prev_P == P_list[ipair]) ? 0 : 1;
		
		write_P = 0;
		if (ipair == npairs - 1)
		{
			write_MN = 1;
			write_P  = 1;
		} else {
			write_MN = 0;
			if (P_list[ipair] != P_list[ipair + 1]) write_P  = 1;
		}
		prev_P = P_list[ipair];
		
		int dimP = TinySCF->shell_bf_num[P_list[ipair]];
		int dimQ = TinySCF->shell_bf_num[Q_list[ipair]];
		int is_1111 = dimM * dimN * dimP * dimQ;
		
		if (is_1111 == 1)    Accum_Fock_1111  (ACCUM_FOCK_PARAM);
		else if (dimQ == 1)  Accum_Fock_dimQ1 (ACCUM_FOCK_PARAM);
		else if (dimQ == 3)  Accum_Fock_dimQ3 (ACCUM_FOCK_PARAM);
		else if (dimQ == 6)  Accum_Fock_dimQ6 (ACCUM_FOCK_PARAM);
		else if (dimQ == 10) Accum_Fock_dimQ10(ACCUM_FOCK_PARAM);
		else if (dimQ == 15) Accum_Fock_dimQ15(ACCUM_FOCK_PARAM);
		else Accum_Fock(ACCUM_FOCK_PARAM);
	}
}

// F = H_core + (J + J^T) / 2 + (K + K^T) / 2;
void TinySCF_HJKmat_to_Fmat(double *Hcore_mat, double *J_mat, double *K_mat, double *F_mat, int nbf)
{
	#pragma omp for
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
	Simint_t simint    = TinySCF->simint;
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
		CMS_Simint_createThreadMultishellpair(&thread_multi_shellpair);
		
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
				int ket_id = CMS_Simint_getShellpairAMIndex(simint, P, Q);
				KetShellpairList_t target_shellpair_list = &thread_ksp_lists->ket_shellpair_lists[ket_id];
				add_shellpair_to_KetShellPairList(target_shellpair_list, P, Q);
				
				// If the ket-side shell pair list we just used is full, handle it
				if (target_shellpair_list->npairs == MAX_LIST_SIZE)
				{
					double *thread_batch_eris;
					int thread_nints;
					
					// Compute batched ERIs
					CMS_computeShellQuartetBatch_Simint(
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
						if (tid == 0) CMS_Simint_addupdateFtimer(simint, et - st);
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
					CMS_computeShellQuartetBatch_Simint(
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
						if (tid == 0) CMS_Simint_addupdateFtimer(simint, et - st);
					}
					
					// Reset the computed ket-side shell pair list
					target_shellpair_list->npairs = 0;
				}
			}
		}
		
		// Free batch ERI auxiliary data structures
		CMS_Simint_freeThreadMultishellpair(&thread_multi_shellpair);
		free_ThreadKetShellpairLists(thread_ksp_lists);
		
		#pragma omp barrier
		TinySCF_HJKmat_to_Fmat(Hcore_mat, J_mat, K_mat, F_mat, num_bas_func);
	}
}
