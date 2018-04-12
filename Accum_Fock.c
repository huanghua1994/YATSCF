#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "TinySCF.h"
#include "Accum_Fock.h"

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

static inline void atomic_add_f64(volatile double *global_value, double addend)
{
	uint64_t expected_value, new_value;
	do {
		double old_value = *global_value;
		expected_value = _castf64_u64(old_value);
		new_value = _castf64_u64(old_value + addend);
	} while (!__sync_bool_compare_and_swap((volatile uint64_t*)global_value, expected_value, new_value));
}

static inline void atomic_update_global_block(
	double *global_block, double *local_block,
	int nrows, int ncols, int ldg
)
{
	for (int irow = 0; irow < nrows; irow++)
	{
		int base1 = irow * ncols;
		int base2 = irow * ldg;
		for (int icol = 0; icol < ncols; icol++)
		{
			double addend = local_block[base1 + icol];
			atomic_add_f64(&global_block[base2 + icol], addend);
		}
	}
}

void Accum_Fock(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = TinySCF->shell_bf_num[Q];
	int idxM = TinySCF->shell_bf_sind[M];
	int idxN = TinySCF->shell_bf_sind[N];
	int idxP = TinySCF->shell_bf_sind[P];
	int idxQ = TinySCF->shell_bf_sind[Q];
	
	// Set global matrix pointers
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
	
	// Set buffer pointer
	double *thread_buf = TinySCF->Accum_Fock_buf + tid * TinySCF->max_buf_size;
	int required_buf_size = (dimP + dimN + dimM) * dimQ + (dimP + dimN + dimM) * dimQ;
	assert(required_buf_size <= TinySCF->max_buf_size);
	double *write_buf = thread_buf;
	double *J_MN_buf = write_buf;  write_buf += dimM * dimN;
	double *K_MP_buf = write_buf;  write_buf += dimM * dimP;
	double *K_NP_buf = write_buf;  write_buf += dimN * dimP;
	double *J_PQ_buf = write_buf;  write_buf += dimP * dimQ;
	double *K_NQ_buf = write_buf;  write_buf += dimN * dimQ;
	double *K_MQ_buf = write_buf;  write_buf += dimM * dimQ;

	// Reset result buffer
	if (load_MN) memset(J_MN_buf, 0, sizeof(double) * dimM * dimN);
	if (load_P)  memset(K_MP_buf, 0, sizeof(double) * dimP * (dimM + dimN));
	memset(J_PQ_buf, 0, sizeof(double) * dimQ * (dimM + dimN + dimP));
	
	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	for (int iM = 0; iM < dimM; iM++) 
	{
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iM_nbf  = iM * nbf;
			int iN_nbf  = iN * nbf;
			int iM_dimQ = iM * dimQ;
			int iN_dimQ = iN * dimQ;
			double coef1_D_MN = coef[1] * D_MN[iM * nbf + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_nbf  = iP * nbf;
				int iP_dimQ = iP * dimQ;
				int Ibase = dimQ * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN * nbf + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM * nbf + iP];
				double k_MP = 0.0, k_NP = 0.0;
				// dimQ is small, vectorizing short loops may hurt performance 
				// since it needs horizon reduction after the loop
				for (int iQ = 0; iQ < dimQ; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_nbf + iQ] * I;
					k_MP -= D_NQ[iN_nbf + iQ] * I;
					k_NP -= D_MQ[iM_nbf + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM * dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN * dimP + iP] += coef[3] * k_NP;
			} // for (int iM = 0; iM < dimM; iM++) 
			J_MN_buf[iM * dimN + iN] += coef[0] * j_MN;
		} // for (int iQ = 0; iQ < dimQ; iQ++) 
	} // for (int iN = 0; iN < dimN; iN++)
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_global_block(J_MN, J_MN_buf, dimM, dimN, nbf);
	
	if (write_P)
	{
		atomic_update_global_block(K_MP, K_MP_buf, dimM, dimP, nbf);
		atomic_update_global_block(K_NP, K_NP_buf, dimN, dimP, nbf);
	}
	
	atomic_update_global_block(J_PQ, J_PQ_buf, dimP, dimQ, nbf);
	atomic_update_global_block(K_MQ, K_MQ_buf, dimM, dimQ, nbf);
	atomic_update_global_block(K_NQ, K_NQ_buf, dimN, dimQ, nbf);
}
