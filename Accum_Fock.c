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

static inline void atomic_update_vector(double *dst, double *src, int length)
{
	for (int i = 0; i < length; i++)
		atomic_add_f64(&dst[i], src[i]);
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
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
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
		int iM_dimP = iM * dimP;
		int iM_dimN = iM * dimN;
		int iM_dimQ = iM * dimQ;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iN_dimP = iN * dimP;
			int iN_dimQ = iN * dimQ;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_dimQ = iP * dimQ;
				int Ibase = dimQ * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];
				double k_MP = 0.0, k_NP = 0.0;
				// dimQ is small, vectorizing short loops may hurt performance 
				// since it needs horizon reduction after the loop
				for (int iQ = 0; iQ < dimQ; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_dimQ + iQ] * I;
					k_MP -= D_NQ[iN_dimQ + iQ] * I;
					k_NP -= D_MQ[iM_dimQ + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM_dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN_dimP + iP] += coef[3] * k_NP;
			}  // for (int iP = 0; iP < dimP; iP++) 
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		} // for (int iN = 0; iN < dimN; iN++) 
	} // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * dimQ);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * dimQ);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * dimQ);
}

// ----- Specialized implementations of Accum_Fock with different dimQ -----
// ----- We don't have function template in C, so we have to copy them -----

void Accum_Fock_dimQ1(
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
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
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
		int iM_dimN = iM * dimN;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iM_dimP = iM * dimP;
			int iN_dimP = iN * dimP;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0, k_MQ = 0.0, k_NQ = 0.0;
			int Ibase = dimP * (iN + dimN * iM);
			for (int iP = 0; iP < dimP; iP++) 
			{
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];

				double I = ERI[Ibase + iP];
				
				j_MN += D_PQ[iP] * I;
				k_MQ += ncoef4_D_NP * I;
				k_NQ += ncoef5_D_MP * I;

				J_PQ_buf[iP] += coef1_D_MN * I;
				K_MP_buf[iM_dimP + iP] -= coef[2] * D_NQ[iN] * I;
				K_NP_buf[iN_dimP + iP] -= coef[3] * D_MQ[iM] * I;
			}  // for (int iP = 0; iP < dimP; iP++)  
			K_MQ_buf[iM] += k_MQ;
			K_NQ_buf[iN] += k_NQ;
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		}  // for (int iN = 0; iN < dimN; iN++) 
	}  // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * dimQ);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * dimQ);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * dimQ);
}

void Accum_Fock_dimQ3(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = 3; // TinySCF->shell_bf_num[Q];
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	// Set buffer pointer
	double *thread_buf = TinySCF->Accum_Fock_buf + tid * TinySCF->max_buf_size;
	int required_buf_size = (dimP + dimN + dimM) * 3 + (dimP + dimN + dimM) * 3;
	assert(required_buf_size <= TinySCF->max_buf_size);
	double *write_buf = thread_buf;
	double *J_MN_buf = write_buf;  write_buf += dimM * dimN;
	double *K_MP_buf = write_buf;  write_buf += dimM * dimP;
	double *K_NP_buf = write_buf;  write_buf += dimN * dimP;
	double *J_PQ_buf = write_buf;  write_buf += dimP * 3;
	double *K_NQ_buf = write_buf;  write_buf += dimN * 3;
	double *K_MQ_buf = write_buf;  write_buf += dimM * 3;

	// Reset result buffer
	if (load_MN) memset(J_MN_buf, 0, sizeof(double) * dimM * dimN);
	if (load_P)  memset(K_MP_buf, 0, sizeof(double) * dimP * (dimM + dimN));
	memset(J_PQ_buf, 0, sizeof(double) * 3 * (dimM + dimN + dimP));
	
	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	for (int iM = 0; iM < dimM; iM++) 
	{
		int iM_dimP = iM * dimP;
		int iM_dimN = iM * dimN;
		int iM_dimQ = iM * 3;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iN_dimP = iN * dimP;
			int iN_dimQ = iN * 3;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_dimQ = iP * 3;
				int Ibase = 3 * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];
				double k_MP = 0.0, k_NP = 0.0;
				
				
				for (int iQ = 0; iQ < 3; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_dimQ + iQ] * I;
					k_MP -= D_NQ[iN_dimQ + iQ] * I;
					k_NP -= D_MQ[iM_dimQ + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM_dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN_dimP + iP] += coef[3] * k_NP;
			}  // for (int iP = 0; iP < dimP; iP++) 
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		}  // for (int iN = 0; iN < dimN; iN++) 
	}  // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * 3);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * 3);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * 3);
}

void Accum_Fock_dimQ6(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = 6; // TinySCF->shell_bf_num[Q];
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	// Set buffer pointer
	double *thread_buf = TinySCF->Accum_Fock_buf + tid * TinySCF->max_buf_size;
	int required_buf_size = (dimP + dimN + dimM) * 6 + (dimP + dimN + dimM) * 6;
	assert(required_buf_size <= TinySCF->max_buf_size);
	double *write_buf = thread_buf;
	double *J_MN_buf = write_buf;  write_buf += dimM * dimN;
	double *K_MP_buf = write_buf;  write_buf += dimM * dimP;
	double *K_NP_buf = write_buf;  write_buf += dimN * dimP;
	double *J_PQ_buf = write_buf;  write_buf += dimP * 6;
	double *K_NQ_buf = write_buf;  write_buf += dimN * 6;
	double *K_MQ_buf = write_buf;  write_buf += dimM * 6;

	// Reset result buffer
	if (load_MN) memset(J_MN_buf, 0, sizeof(double) * dimM * dimN);
	if (load_P)  memset(K_MP_buf, 0, sizeof(double) * dimP * (dimM + dimN));
	memset(J_PQ_buf, 0, sizeof(double) * 6 * (dimM + dimN + dimP));
	
	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	for (int iM = 0; iM < dimM; iM++) 
	{
		int iM_dimP = iM * dimP;
		int iM_dimN = iM * dimN;
		int iM_dimQ = iM * 6;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iN_dimP = iN * dimP;
			int iN_dimQ = iN * 6;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_dimQ = iP * 6;
				int Ibase = 6 * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];
				double k_MP = 0.0, k_NP = 0.0;
				
				for (int iQ = 0; iQ < 6; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_dimQ + iQ] * I;
					k_MP -= D_NQ[iN_dimQ + iQ] * I;
					k_NP -= D_MQ[iM_dimQ + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM_dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN_dimP + iP] += coef[3] * k_NP;
			}  // for (int iP = 0; iP < dimP; iP++) 
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		} // for (int iN = 0; iN < dimN; iN++) 
	} // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * 6);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * 6);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * 6);
}

void Accum_Fock_dimQ10(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = 10; // TinySCF->shell_bf_num[Q];
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	// Set buffer pointer
	double *thread_buf = TinySCF->Accum_Fock_buf + tid * TinySCF->max_buf_size;
	int required_buf_size = (dimP + dimN + dimM) * 10 + (dimP + dimN + dimM) * 10;
	assert(required_buf_size <= TinySCF->max_buf_size);
	double *write_buf = thread_buf;
	double *J_MN_buf = write_buf;  write_buf += dimM * dimN;
	double *K_MP_buf = write_buf;  write_buf += dimM * dimP;
	double *K_NP_buf = write_buf;  write_buf += dimN * dimP;
	double *J_PQ_buf = write_buf;  write_buf += dimP * 10;
	double *K_NQ_buf = write_buf;  write_buf += dimN * 10;
	double *K_MQ_buf = write_buf;  write_buf += dimM * 10;

	// Reset result buffer
	if (load_MN) memset(J_MN_buf, 0, sizeof(double) * dimM * dimN);
	if (load_P)  memset(K_MP_buf, 0, sizeof(double) * dimP * (dimM + dimN));
	memset(J_PQ_buf, 0, sizeof(double) * 10 * (dimM + dimN + dimP));
	
	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	for (int iM = 0; iM < dimM; iM++) 
	{
		int iM_dimP = iM * dimP;
		int iM_dimN = iM * dimN;
		int iM_dimQ = iM * 10;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iN_dimP = iN * dimP;
			int iN_dimQ = iN * 10;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_dimQ = iP * 10;
				int Ibase = 10 * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];
				double k_MP = 0.0, k_NP = 0.0;
				
				for (int iQ = 0; iQ < 10; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_dimQ + iQ] * I;
					k_MP -= D_NQ[iN_dimQ + iQ] * I;
					k_NP -= D_MQ[iM_dimQ + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM_dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN_dimP + iP] += coef[3] * k_NP;
			}  // for (int iP = 0; iP < dimP; iP++) 
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		} // for (int iN = 0; iN < dimN; iN++) 
	} // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * dimQ);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * dimQ);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * dimQ);
}

void Accum_Fock_dimQ15(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nbf  = TinySCF->nbasfuncs;
	int dimM = TinySCF->shell_bf_num[M];
	int dimN = TinySCF->shell_bf_num[N];
	int dimP = TinySCF->shell_bf_num[P];
	int dimQ = 15; // TinySCF->shell_bf_num[Q];
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	// Set buffer pointer
	double *thread_buf = TinySCF->Accum_Fock_buf + tid * TinySCF->max_buf_size;
	int required_buf_size = (dimP + dimN + dimM) * 15 + (dimP + dimN + dimM) * 15;
	assert(required_buf_size <= TinySCF->max_buf_size);
	double *write_buf = thread_buf;
	double *J_MN_buf = write_buf;  write_buf += dimM * dimN;
	double *K_MP_buf = write_buf;  write_buf += dimM * dimP;
	double *K_NP_buf = write_buf;  write_buf += dimN * dimP;
	double *J_PQ_buf = write_buf;  write_buf += dimP * 15;
	double *K_NQ_buf = write_buf;  write_buf += dimN * 15;
	double *K_MQ_buf = write_buf;  write_buf += dimM * 15;

	// Reset result buffer
	if (load_MN) memset(J_MN_buf, 0, sizeof(double) * dimM * dimN);
	if (load_P)  memset(K_MP_buf, 0, sizeof(double) * dimP * (dimM + dimN));
	memset(J_PQ_buf, 0, sizeof(double) * 15 * (dimM + dimN + dimP));
	
	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	for (int iM = 0; iM < dimM; iM++) 
	{
		int iM_dimP = iM * dimP;
		int iM_dimN = iM * dimN;
		int iM_dimQ = iM * 15;
		for (int iN = 0; iN < dimN; iN++) 
		{
			int iN_dimP = iN * dimP;
			int iN_dimQ = iN * 15;
			double coef1_D_MN = coef[1] * D_MN[iM_dimN + iN];
			double j_MN = 0.0;
			for (int iP = 0; iP < dimP; iP++) 
			{
				int iP_dimQ = iP * 15;
				int Ibase = dimQ * (iP + dimP * (iN + dimN * iM));
				double ncoef4_D_NP = -coef[4] * D_NP[iN_dimP + iP];
				double ncoef5_D_MP = -coef[5] * D_MP[iM_dimP + iP];
				double k_MP = 0.0, k_NP = 0.0;
				
				#pragma simd
				for (int iQ = 0; iQ < 15; iQ++) 
				{
					double I = ERI[Ibase + iQ];
					
					j_MN += D_PQ[iP_dimQ + iQ] * I;
					k_MP -= D_NQ[iN_dimQ + iQ] * I;
					k_NP -= D_MQ[iM_dimQ + iQ] * I;

					J_PQ_buf[iP_dimQ + iQ] +=  coef1_D_MN * I;
					K_MQ_buf[iM_dimQ + iQ] += ncoef4_D_NP * I;
					K_NQ_buf[iN_dimQ + iQ] += ncoef5_D_MP * I;
				}
				K_MP_buf[iM_dimP + iP] += coef[2] * k_MP;
				K_NP_buf[iN_dimP + iP] += coef[3] * k_NP;
			}  // for (int iP = 0; iP < dimP; iP++) 
			J_MN_buf[iM_dimN + iN] += coef[0] * j_MN;
		} // for (int iN = 0; iN < dimN; iN++) 
	} // for (int iM = 0; iM < dimM; iM++) 
	
	// Update to global array using atomic_add_f64()
	if (write_MN) atomic_update_vector(J_MN, J_MN_buf, dimM * dimN);
	
	if (write_P)
	{
		atomic_update_vector(K_MP, K_MP_buf, dimM * dimP);
		atomic_update_vector(K_NP, K_NP_buf, dimN * dimP);
	}
	
	atomic_update_vector(J_PQ, J_PQ_buf, dimP * 15);
	atomic_update_vector(K_MQ, K_MQ_buf, dimM * 15);
	atomic_update_vector(K_NQ, K_NQ_buf, dimN * 15);
}

void Accum_Fock_1111(
	TinySCF_t TinySCF, int tid, int M, int N, int P, int Q, double *ERI,
	int load_MN, int load_P, int write_MN, int write_P
)
{
	// Set matrix size info
	int nshells = TinySCF->nshells;
	
	// Set global matrix pointers
	double *J_MN = TinySCF->J_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *J_PQ = TinySCF->J_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *K_MP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *K_NP = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *K_MQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *K_NQ = TinySCF->K_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];
	
	double *D_MN = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + N];
	double *D_PQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[P * nshells + Q];
	double *D_MP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + P];
	double *D_NP = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + P];
	double *D_MQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[M * nshells + Q];
	double *D_NQ = TinySCF->D_mat_block + TinySCF->mat_block_ptr[N * nshells + Q];

	// Get uniqueness ERI symmetric 
	double coef[7];
	unique_integral_coef(M, N, P, Q, coef);
	
	double I = ERI[0];
	
	double vMN =  coef[0] * D_PQ[0] * I;
	double vPQ =  coef[1] * D_MN[0] * I;
	double vMP = -coef[2] * D_NQ[0] * I;
	double vNP = -coef[3] * D_MQ[0] * I;
	double vMQ = -coef[4] * D_NP[0] * I;
	double vNQ = -coef[5] * D_MP[0] * I;
	
	atomic_add_f64(&J_MN[0], vMN);
	atomic_add_f64(&J_PQ[0], vPQ);
	atomic_add_f64(&K_MP[0], vMP);
	atomic_add_f64(&K_NP[0], vNP);
	atomic_add_f64(&K_MQ[0], vMQ);
	atomic_add_f64(&K_NQ[0], vNQ);
}