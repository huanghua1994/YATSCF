#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "utils.h"
#include "TinySCF.h"
#include "build_Fock.h"

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
	double *J_mat, double *K_mat, double *D_mat, double *ERI, int tid, 
	int nbf, int M, int N, int P, int Q, int *shell_bf_num, int *shell_bf_sind
)
{
	int dimM = shell_bf_num[M];
	int dimN = shell_bf_num[N];
	int dimP = shell_bf_num[P];
	int dimQ = shell_bf_num[Q];
	int idxM = shell_bf_sind[M];
	int idxN = shell_bf_sind[N];
	int idxP = shell_bf_sind[P];
	int idxQ = shell_bf_sind[Q];
	
	double *J_MN = J_mat + idxM * nbf + idxN;
	double *J_PQ = J_mat + idxP * nbf + idxQ;
	double *K_MP = K_mat + idxM * nbf + idxP;
	double *K_NP = K_mat + idxN * nbf + idxP;
	double *K_MQ = K_mat + idxM * nbf + idxQ;
	double *K_NQ = K_mat + idxN * nbf + idxQ;
	
	double *D_MN = D_mat + idxM * nbf + idxN;
	double *D_PQ = D_mat + idxP * nbf + idxQ;
	double *D_MP = D_mat + idxM * nbf + idxP;
	double *D_NP = D_mat + idxN * nbf + idxP;
	double *D_MQ = D_mat + idxM * nbf + idxQ;
	double *D_NQ = D_mat + idxN * nbf + idxQ;
	
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
				
				int nints;
				double *integrals;
				CInt_computeShellQuartet_SIMINT(simint, tid, M, N, P, Q, &integrals, &nints);
				double st = get_wtime_sec();
				Accum_Fock(
					J_mat, K_mat, D_mat, integrals, tid, num_bas_func, 
					M, N, P, Q, shell_bf_num, shell_bf_sind
				);
				double et = get_wtime_sec();
				if (tid == 0) CInt_SIMINT_addupdateFtimer(simint, et - st);
			}
		}
	}
	
	TinySCF_HJKmat_to_Fmat(Hcore_mat, J_mat, K_mat, F_mat, num_bas_func);
}
