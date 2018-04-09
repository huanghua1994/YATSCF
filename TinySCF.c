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
	TinySCF->prim_screen_tol2  = 1e-14 * 1e-14;
	TinySCF->shell_screen_tol2 = 1e-11 * 1e-11;
	
	// Allocate memory for matrices
	TinySCF->sp_scr_vals = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->nshellpairs);
	TinySCF->H_core      = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->Ovlp_mat    = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->Fock_mat    = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size);
	TinySCF->prev_Fock   = (double*) ALIGN64B_MALLOC(DBL_SIZE * TinySCF->mat_size * MAX_DIIS);
	assert(TinySCF->sp_scr_vals   != NULL);
	assert(TinySCF->H_core        != NULL);
	assert(TinySCF->Ovlp_mat      != NULL);
	assert(TinySCF->Fock_mat      != NULL);
	assert(TinySCF->prev_Fock     != NULL);
	
	// Allocate memory for all shells' basis function info
	TinySCF->shell_bf_sind = (int*) ALIGN64B_MALLOC(INT_SIZE * (TinySCF->nshells + 1));
	TinySCF->shell_bf_num  = (int*) ALIGN64B_MALLOC(INT_SIZE * TinySCF->nshells);
	assert(TinySCF->shell_bf_sind != NULL);
	assert(TinySCF->shell_bf_num  != NULL);
	for (int i = 0; i < TinySCF->nshells; i++)
	{
		TinySCF->shell_bf_sind[i] = CInt_getFuncStartInd(TinySCF->basis, i);
		TinySCF->shell_bf_num[i]  = CInt_getShellDim    (TinySCF->basis, i);
		//printf("[DEBUG] Shell %d: %d basis functions, basis function start index = %d\n", 
		//		i, TinySCF->shell_bf_num[i], TinySCF->shell_bf_sind[i]);
	}
	
	TinySCF->mem_size      = DBL_SIZE * (TinySCF->nshellpairs + (MAX_DIIS + 3) * TinySCF->mat_size) 
	                       + INT_SIZE * (2 * TinySCF->nshells + 1);
	printf("TinySCF initialization over, memory usage: %.2lf MB\n", (double) TinySCF->mem_size / 1048576.0);
}

void free_TinySCF(TinySCF_t TinySCF)
{
	ALIGN64B_FREE(TinySCF->sp_scr_vals);
	ALIGN64B_FREE(TinySCF->H_core);
	ALIGN64B_FREE(TinySCF->Ovlp_mat);
	ALIGN64B_FREE(TinySCF->Fock_mat);
	ALIGN64B_FREE(TinySCF->prev_Fock);
	
	ALIGN64B_FREE(TinySCF->shell_bf_sind);
	ALIGN64B_FREE(TinySCF->shell_bf_num);
	
	free(TinySCF);
}