#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mkl.h>

#include "CInt.h"
#include "utils.h"
#include "TinySCF.h"
#include "build_Fock.h"
#include "build_density.h"
#include "DIIS.h"

// This file only contains functions that initializing TinySCF engine, precomputing reusable 
// matrices and arrays and destroying TinySCF engine. Most time consuming functions are in
// build_Fock.c, build_density.c and DIIS.c

void init_TinySCF(TinySCF_t TinySCF, char *bas_fname, char *xyz_fname, const int niters)
{
	assert(TinySCF != NULL);
	
	double st = get_wtime_sec();
	
	TinySCF->nthreads = omp_get_max_threads();
	
	// Reset statistic info
	TinySCF->mem_size       = 0.0;
	TinySCF->init_time      = 0.0;
	TinySCF->S_Hcore_time   = 0.0;
	TinySCF->shell_scr_time = 0.0;
	
	// Load basis set and molecule from input
	CInt_createBasisSet(&(TinySCF->basis));
	CInt_loadBasisSet(TinySCF->basis, bas_fname, xyz_fname);
	
	// Calculate problem info
	TinySCF->natoms      = CInt_getNumAtoms   (TinySCF->basis);
	TinySCF->nshells     = CInt_getNumShells  (TinySCF->basis);
	TinySCF->nbasfuncs   = CInt_getNumFuncs   (TinySCF->basis);
	TinySCF->charge      = CInt_getTotalCharge(TinySCF->basis);
	TinySCF->n_occ       = CInt_getNumOccOrb  (TinySCF->basis);
	TinySCF->niters      = niters;
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
	TinySCF->shell_scrtol2 = 1e-11 * 1e-11;
	
	// Set SCF iteration termination criteria 
	TinySCF->energy_delta_tol = 1e-10;
	
	// Allocate memory for matrices
	int mat_mem_size = DBL_SIZE * TinySCF->mat_size;
	TinySCF->Hcore_mat   = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->S_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->F_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->D_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->J_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->K_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->X_mat       = (double*) ALIGN64B_MALLOC(mat_mem_size);
	TinySCF->prev_F_mat  = (double*) ALIGN64B_MALLOC(mat_mem_size * MAX_DIIS);
	assert(TinySCF->Hcore_mat     != NULL);
	assert(TinySCF->S_mat         != NULL);
	assert(TinySCF->F_mat         != NULL);
	assert(TinySCF->D_mat         != NULL);
	assert(TinySCF->J_mat         != NULL);
	assert(TinySCF->K_mat         != NULL);
	assert(TinySCF->X_mat         != NULL);
	assert(TinySCF->prev_F_mat    != NULL);
	TinySCF->mem_size += (MAX_DIIS + 7) * mat_mem_size;

	// Allocate memory for temporary arrays in building density matrix
	TinySCF->Cocc_mat    = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->n_occ * TinySCF->nbasfuncs);
	TinySCF->eigval      = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->nbasfuncs);
	TinySCF->ev_idx      = (int*)    ALIGN64B_MALLOC(INT_SIZE * TinySCF->nbasfuncs);
	assert(TinySCF->Cocc_mat      != NULL);
	assert(TinySCF->eigval        != NULL);
	assert(TinySCF->eigval        != NULL);
	TinySCF->mem_size += DBL_SIZE * (TinySCF->n_occ + 1) * TinySCF->nbasfuncs;
	TinySCF->mem_size += INT_SIZE * TinySCF->nbasfuncs;
	
	// Allocate memory for all shells' basis function info
	TinySCF->num_uniq_sp   = (TinySCF->nshells + 1) * TinySCF->nshells / 2;
	TinySCF->sp_scrval     = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->nshellpairs);
	TinySCF->shell_bf_sind = (int*)    ALIGN64B_MALLOC(INT_SIZE * (TinySCF->nshells + 1));
	TinySCF->shell_bf_num  = (int*)    ALIGN64B_MALLOC(INT_SIZE * TinySCF->nshells);
	TinySCF->uniq_sp_lid   = (int*)    ALIGN64B_MALLOC(INT_SIZE * TinySCF->num_uniq_sp);
	TinySCF->uniq_sp_rid   = (int*)    ALIGN64B_MALLOC(INT_SIZE * TinySCF->num_uniq_sp);
	assert(TinySCF->sp_scrval     != NULL);
	assert(TinySCF->shell_bf_sind != NULL);
	assert(TinySCF->shell_bf_num  != NULL);
	assert(TinySCF->uniq_sp_lid   != NULL);
	assert(TinySCF->uniq_sp_rid   != NULL);
	TinySCF->mem_size += DBL_SIZE * TinySCF->nshellpairs;
	TinySCF->mem_size += INT_SIZE * (2 * TinySCF->nshells + 1);
	TinySCF->mem_size += INT_SIZE * 2 * TinySCF->num_uniq_sp;
	for (int i = 0; i < TinySCF->nshells; i++)
	{
		TinySCF->shell_bf_sind[i] = CInt_getFuncStartInd(TinySCF->basis, i);
		TinySCF->shell_bf_num[i]  = CInt_getShellDim    (TinySCF->basis, i);
	}
	TinySCF->shell_bf_sind[TinySCF->nshells] = TinySCF->nbasfuncs;
	
	// Initialize Simint object
	CInt_createSIMINT(TinySCF->basis, &(TinySCF->simint), TinySCF->nthreads);
	
	double et = get_wtime_sec();
	TinySCF->init_time = et - st;
	
	// Print memory usage
	printf("TinySCF initialization over, memory usage: %.2lf MB, ", TinySCF->mem_size / 1048576.0);
	printf("time elapsed = %.2lf (s)\n", TinySCF->init_time);
}


void free_TinySCF(TinySCF_t TinySCF)
{
	assert(TinySCF != NULL);
	
	ALIGN64B_FREE(TinySCF->Hcore_mat);
	ALIGN64B_FREE(TinySCF->S_mat);
	ALIGN64B_FREE(TinySCF->F_mat);
	ALIGN64B_FREE(TinySCF->D_mat);
	ALIGN64B_FREE(TinySCF->J_mat);
	ALIGN64B_FREE(TinySCF->K_mat);
	ALIGN64B_FREE(TinySCF->X_mat);
	ALIGN64B_FREE(TinySCF->prev_F_mat);
	
	ALIGN64B_FREE(TinySCF->Cocc_mat);
	ALIGN64B_FREE(TinySCF->eigval);
	ALIGN64B_FREE(TinySCF->ev_idx);
	
	ALIGN64B_FREE(TinySCF->sp_scrval);
	ALIGN64B_FREE(TinySCF->shell_bf_sind);
	ALIGN64B_FREE(TinySCF->shell_bf_num);
	ALIGN64B_FREE(TinySCF->uniq_sp_lid);
	ALIGN64B_FREE(TinySCF->uniq_sp_rid);
	
	CInt_destroyBasisSet(TinySCF->basis);
	CInt_destroySIMINT(TinySCF->simint, 1);
	
	free(TinySCF);
}

void TinySCF_compute_Hcore_Ovlp_mat(TinySCF_t TinySCF)
{
	assert(TinySCF != NULL);
	
	double st = get_wtime_sec();
	
	// Compute core Hamiltonian and overlap matrix
	memset(TinySCF->Hcore_mat, 0, DBL_SIZE * TinySCF->mat_size);
	memset(TinySCF->S_mat,     0, DBL_SIZE * TinySCF->mat_size);
	for (int M = 0; M < TinySCF->nshells; M++)
	{
		for (int N = 0; N < TinySCF->nshells; N++)
		{
			int nints, tid = 0;
			double *integrals;
			
			int mat_topleft_offset = TinySCF->shell_bf_sind[M] * TinySCF->nbasfuncs + TinySCF->shell_bf_sind[N];
			double *S_mat_ptr      = TinySCF->S_mat     + mat_topleft_offset;
			double *Hcore_mat_ptr  = TinySCF->Hcore_mat + mat_topleft_offset;
			
			int nrows = TinySCF->shell_bf_num[M];
			int ncols = TinySCF->shell_bf_num[N];
			
			// Compute the contribution of current shell pair to core Hamiltonian matrix
			CInt_computePairOvl_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			if (nints > 0) copy_matrix_block(S_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
			
			// Compute the contribution of current shell pair to overlap matrix
			CInt_computePairCoreH_SIMINT(TinySCF->basis, TinySCF->simint, tid, M, N, &integrals, &nints);
			if (nints > 0) copy_matrix_block(Hcore_mat_ptr, TinySCF->nbasfuncs, integrals, ncols, nrows, ncols);
		}
	}
	
	// Construct basis transformation 
	int N = TinySCF->nbasfuncs;
	double *U_mat  = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	double *U_mat0 = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	double *eigval = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->nbasfuncs);
	assert(U_mat != NULL && eigval != NULL);
	memcpy(U_mat, TinySCF->S_mat, DBL_SIZE * TinySCF->mat_size);
	// [U, D] = eig(S);
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, U_mat, N, eigval); 
	// X = U * D^{-1/2} * U'
	memcpy(U_mat0, U_mat, DBL_SIZE * TinySCF->mat_size);
	for (int i = 0; i < N; i++) 
		eigval[i] = 1.0 / sqrt(eigval[i]);
	for (int i = 0; i < N; i++)
	{
		#pragma simd
		for (int j = 0; j < N; j++)
			U_mat0[i * N + j] *= eigval[j];
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, U_mat0, N, U_mat, N, 0.0, TinySCF->X_mat, N);
	
	ALIGN64B_FREE(U_mat);
	ALIGN64B_FREE(U_mat0);
	ALIGN64B_FREE(eigval);
	
	double et = get_wtime_sec();
	TinySCF->S_Hcore_time = et - st;
	
	// Print runtime
	printf("TinySCF precompute Hcore, S, and X matrices over,   ");
	printf("time elapsed = %.2lf (s)\n", TinySCF->S_Hcore_time);
}

void TinySCF_compute_sq_Schwarz_scrvals(TinySCF_t TinySCF)
{
	assert(TinySCF != NULL);
	
	double st = get_wtime_sec();
	
	// Compute screening values using Schwarz inequality
	double global_max_scrval = 0.0;
	for (int M = 0; M < TinySCF->nshells; M++)
	{
		int dimM = TinySCF->shell_bf_num[M];
		for (int N = 0; N < TinySCF->nshells; N++)
		{
			int dimN = TinySCF->shell_bf_num[N];
			
			int nints;
			double *integrals;
			CInt_computeShellQuartet_SIMINT(TinySCF->simint, 0, M, N, M, N, &integrals, &nints);
			
			double maxval = 0.0;
			if (nints > 0)
			{
				// Loop over all ERIs in a shell quartet and find the max value
				for (int iM = 0; iM < dimM; iM++)
					for (int iN = 0; iN < dimN; iN++)
					{
						int index = iN * (dimM * dimN * dimM + dimM) + iM * (dimN * dimM + 1); // Simint layout
						double val = fabs(integrals[index]);
						if (val > maxval) maxval = val;
					}
			}
			TinySCF->sp_scrval[M * TinySCF->nshells + N] = maxval;
			if (maxval > global_max_scrval) global_max_scrval = maxval;
		}
	}
	TinySCF->max_scrval = global_max_scrval;
	
	// Recreate Simint object to reset its statistic info
	CInt_destroySIMINT(TinySCF->simint, 0);
	CInt_createSIMINT(TinySCF->basis, &(TinySCF->simint), TinySCF->nthreads);
	
	// Generate unique shell pairs that survive Schwarz screening
	// eta is the threshold for screening a shell pair
	double eta = TinySCF->shell_scrtol2 / TinySCF->max_scrval;
	int nnz = 0;
	for (int M = 0; M < TinySCF->nshells; M++)
	{
		for (int N = 0; N < TinySCF->nshells; N++)
		{
			double sp_scrval = TinySCF->sp_scrval[M * TinySCF->nshells + N];
			// if sp_scrval * max_scrval < shell_scrtol2, for any given shell pair
			// (P,Q), (MN|PQ) is always < shell_scrtol2 and will be screened
			if (sp_scrval > eta)  
			{
				// Symmetric uniqueness check, from GTFock, not so sure about the mechanism...
				if ((M > N) && ((M + N) % 2 == 1)) continue;
				if ((M < N) && ((M + N) % 2 == 0)) continue;
				
				TinySCF->uniq_sp_lid[nnz] = M;
				TinySCF->uniq_sp_rid[nnz] = N;
				nnz++;
			}
		}
	}
	TinySCF->num_uniq_sp = nnz;
	
	double et = get_wtime_sec();
	TinySCF->shell_scr_time = et - st;
	
	// Print runtime
	printf("TinySCF precompute shell screening info over,       ");
	printf("time elapsed = %.2lf (s)\n", TinySCF->shell_scr_time);
}

static void TinySCF_get_initial_guess(TinySCF_t TinySCF)
{
	// TODO: use CInt_getInitialGuess() to get initial guess
	memset(TinySCF->D_mat, 0, DBL_SIZE * TinySCF->mat_size);
}

void TinySCF_calc_energy(TinySCF_t TinySCF)
{
	double energy = 0.0;
	#pragma simd 
	for (int i = 0; i < TinySCF->mat_size; i++)
		energy += TinySCF->D_mat[i] * (TinySCF->F_mat[i] + TinySCF->Hcore_mat[i]);
	
	TinySCF->energy = TinySCF->nuc_energy + energy;
}

void TinySCF_do_SCF(TinySCF_t TinySCF)
{
	// Get initial guess for density matrix
	TinySCF_get_initial_guess(TinySCF);
	
	// Calculate nuclear energy
	TinySCF->nuc_energy = CInt_getNucEnergy(TinySCF->basis);
	
	// Start SCF iterations
	printf("TinySCF SCF iteration started...\n");
	TinySCF->iter = 0;
	double prev_energy  = 0;
	double energy_delta = 2018;
	while ((TinySCF->iter < TinySCF->niters) && (energy_delta >= TinySCF->energy_delta_tol))
	{
		printf("--------------- Iteration %d ---------------\n", TinySCF->iter);
		
		double st0, et0, st1, et1;
		st0 = get_wtime_sec();
		
		// Build the Fock matrix
		st1 = get_wtime_sec();
		TinySCF_build_FockMat(TinySCF);
		et1 = get_wtime_sec();
		printf("* Build Fock matrix     : %.2lf (s)\n", et1 - st1);
		
		// DIIS (Pulay mixing)
		st1 = get_wtime_sec();
		TinySCF_DIIS(TinySCF);
		et1 = get_wtime_sec();
		printf("* DIIS procedure        : %.2lf (s)\n", et1 - st1);
		
		// Calculate new system energy
		st1 = get_wtime_sec();
		TinySCF_calc_energy(TinySCF);
		et1 = get_wtime_sec();
		printf("* Calculate energy      : %.2lf (s)\n", et1 - st1);
		
		energy_delta = fabs(TinySCF->energy - prev_energy);
		prev_energy = TinySCF->energy;
		
		// Diagonalize and build the density matrix
		st1 = get_wtime_sec();
		TinySCF_build_DenMat(TinySCF);
		et1 = get_wtime_sec();
		printf("* Build density matrix  : %.2lf (s)\n", et1 - st1);
		
		TinySCF->iter++;
		
		et0 = get_wtime_sec();
		
		printf(
			"* Energy = %.10lf (%.10lf), delta = %.10lf, iteration runtime = %.2lf (s)\n", 
			TinySCF->energy, TinySCF->energy - TinySCF->nuc_energy, energy_delta, et0 - st0
		);
	}
}

