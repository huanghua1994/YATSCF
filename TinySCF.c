#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "CInt.h"
#include "TinySCF.h"

#define ALIGN64B_MALLOC(x) _mm_malloc((x), 64)
#define ALIGN64B_FREE(x)   _mm_free(x)
#define DBL_SIZE           sizeof(double)
#define INT_SIZE           sizeof(int)

void init_TinySCF(TinySCF_t TinySCF, char *bas_fname, char *xyz_fname, const int niters)
{
	assert(TinySCF != NULL);
	
	TinySCF->nthreads = omp_get_max_threads();
	TinySCF->mem_size = 0;
	
	// Load basis set and molecule from input
	BasisSet_t _basis;
	CInt_createBasisSet(&(TinySCF->basis));
	CInt_loadBasisSet(TinySCF->basis, bas_fname, xyz_fname);
	
	// Calculate problem info
	TinySCF->natoms    = CInt_getNumAtoms   (TinySCF->basis);
	TinySCF->nshells   = CInt_getNumShells  (TinySCF->basis);
	TinySCF->nbasfuncs = CInt_getNumFuncs   (TinySCF->basis);
	TinySCF->charge    = CInt_getTotalCharge(TinySCF->basis);
	TinySCF->niters    = niters;
	
	TinySCF->nshellpairs = TinySCF->nshells   * TinySCF->nshells;
	TinySCF->mat_size    = TinySCF->nbasfuncs * TinySCF->nbasfuncs;
	
	printf("Job information:\n");
	printf("    basis set file    = %s\n", bas_fname);
	printf("    molecule  file    = %s\n", xyz_fname);
	printf("    charge            = %d\n", TinySCF->charge);
	printf("    # atoms           = %d\n", TinySCF->natoms);
	printf("    # shells          = %d\n", TinySCF->nshells);
	printf("    # basis functions = %d\n", TinySCF->nbasfuncs);
	
	// Set screening thresholds
	TinySCF->shell_screen_tol2 = 1e-11 * 1e-11;
	
	// Allocate memory for matrices
	TinySCF->sp_scr_vals = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->nshellpairs);
	TinySCF->Hcore_mat   = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->Ovlp_mat    = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->Fock_mat    = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->prev_Fock   = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size * MAX_DIIS);
	assert(TinySCF->sp_scr_vals   != NULL);
	assert(TinySCF->Hcore_mat     != NULL);
	assert(TinySCF->Ovlp_mat      != NULL);
	assert(TinySCF->Fock_mat      != NULL);
	assert(TinySCF->prev_Fock     != NULL);
	TinySCF->mem_size += DBL_SIZE * (TinySCF->nshellpairs + (MAX_DIIS + 3) * TinySCF->mat_size);

	// Allocate memory for all shells' basis function info
	TinySCF->shell_bf_sind = (int*) ALIGN64B_MALLOC(INT_SIZE * (TinySCF->nshells + 1));
	TinySCF->shell_bf_num  = (int*) ALIGN64B_MALLOC(INT_SIZE * TinySCF->nshells);
	assert(TinySCF->shell_bf_sind != NULL);
	assert(TinySCF->shell_bf_num  != NULL);
	TinySCF->mem_size += INT_SIZE * (2 * TinySCF->nshells + 1);
	for (int i = 0; i < TinySCF->nshells; i++)
	{
		TinySCF->shell_bf_sind[i] = CInt_getFuncStartInd(TinySCF->basis, i);
		TinySCF->shell_bf_num[i]  = CInt_getShellDim    (TinySCF->basis, i);
	}
	TinySCF->shell_bf_sind[TinySCF->nshells] = TinySCF->nbasfuncs;
	
	// Initialize Simint object
	CInt_createSIMINT(TinySCF->basis, &(TinySCF->simint), TinySCF->nthreads);
	
	// Print memory usage
	printf("TinySCF initialization over, memory usage: %.2lf MB\n", (double) TinySCF->mem_size / 1048576.0);
}

static void copy_matrix_block(
	double *dst, const int ldd, double *src, const int lds, 
	const int nrows, const int ncols
)
{
	for (int irow = 0; irow < nrows; irow++)
	{
		double *dst_ptr = dst + irow * ldd;
		double *src_ptr = src + irow * lds;
		#pragma simd
		for (int icol = 0; icol < ncols; icol++)
			dst_ptr[icol] = src_ptr[icol];
	}
} 

void TinySCF_compute_Hcore_Ovlp_mat(TinySCF_t TinySCF)
{
	assert(TinySCF != NULL);
	
	memset(TinySCF->Hcore_mat, 0, DBL_SIZE * TinySCF->mat_size);
	memset(TinySCF->Ovlp_mat,  0, DBL_SIZE * TinySCF->mat_size);

	for (int M = 0; M < TinySCF->nshells; M++)
		for (int N = 0; N < TinySCF->nshells; N++)
		{
			int nints, tid = 0;
			double *integrals;
			
			int mat_topleft_offset = TinySCF->shell_bf_sind[M] * TinySCF->nbasfuncs + TinySCF->shell_bf_sind[N];
			double *Ovlp_mat_ptr   = TinySCF->Ovlp_mat  + mat_topleft_offset;
			double *Hcore_mat_ptr  = TinySCF->Hcore_mat + mat_topleft_offset;
			
			int nrows = TinySCF->shell_bf_sind[M + 1] - TinySCF->shell_bf_sind[M];
			int ncols = TinySCF->shell_bf_sind[N + 1] - TinySCF->shell_bf_sind[N];
			
			// Compute the contribution of current shell pair to core Hamiltonian matrix
			//CInt_computePairOvl_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			//if (nints > 0) copy_matrix_block(Ovlp_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
			
			// Compute the contribution of current shell pair to overlap matrix
			CInt_computePairCoreH_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			if (nints > 0) copy_matrix_block(Hcore_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
		}
	
	for (int M = 0; M < TinySCF->nshells; M++)
		for (int N = 0; N < TinySCF->nshells; N++)
		{
			int nints, tid = 0;
			double *integrals;
			
			int mat_topleft_offset = TinySCF->shell_bf_sind[M] * TinySCF->nbasfuncs + TinySCF->shell_bf_sind[N];
			double *Ovlp_mat_ptr   = TinySCF->Ovlp_mat  + mat_topleft_offset;
			double *Hcore_mat_ptr  = TinySCF->Hcore_mat + mat_topleft_offset;
			
			int nrows = TinySCF->shell_bf_sind[M + 1] - TinySCF->shell_bf_sind[M];
			int ncols = TinySCF->shell_bf_sind[N + 1] - TinySCF->shell_bf_sind[N];
			
			// Compute the contribution of current shell pair to core Hamiltonian matrix
			CInt_computePairOvl_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			if (nints > 0) copy_matrix_block(Ovlp_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
			
			// Compute the contribution of current shell pair to overlap matrix
			//CInt_computePairCoreH_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			//if (nints > 0) copy_matrix_block(Hcore_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
		}
	
	// For debug
	// Print overlap matrix
	printf("Overlap:\n");
	for (int i = 0; i < TinySCF->nbasfuncs; i++)
	{
		for (int j = 0; j < TinySCF->nbasfuncs; j++) 
		{
			int idx = i * TinySCF->nbasfuncs + j;
			double x = TinySCF->Ovlp_mat[idx];
			if (x >= 0.0) printf(" ");
			printf("%.7lf\t", x);
		}
		printf("\n");
	}
	printf("\n");
	// Print core Hamiltonian matrix
	printf("core Hamiltonian:\n");
	for (int i = 0; i < TinySCF->nbasfuncs; i++)
	{
		for (int j = 0; j < TinySCF->nbasfuncs; j++) 
		{
			int idx = i * TinySCF->nbasfuncs + j;
			double x = TinySCF->Hcore_mat[idx];
			if (x >= 0.0) printf(" ");
			printf("%.7lf\t", x);
		}
		printf("\n");
	}
	printf("\n");
}

void free_TinySCF(TinySCF_t TinySCF)
{
	assert(TinySCF != NULL);
	
	ALIGN64B_FREE(TinySCF->sp_scr_vals);
	ALIGN64B_FREE(TinySCF->Hcore_mat);
	ALIGN64B_FREE(TinySCF->Ovlp_mat);
	ALIGN64B_FREE(TinySCF->Fock_mat);
	ALIGN64B_FREE(TinySCF->prev_Fock);
	ALIGN64B_FREE(TinySCF->shell_bf_sind);
	ALIGN64B_FREE(TinySCF->shell_bf_num);
	
	free(TinySCF);
}
