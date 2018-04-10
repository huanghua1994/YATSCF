#ifndef _YATSCF_TFOCK_H_
#define _YATSCF_TFOCK_H_

#include <omp.h>
#include "CInt.h"

#define MAX_DIIS 10
#define MIN_DIIS 2

// Tiny SCF engine
struct TinySCF_struct 
{
	int nthreads;
	
	// Problem info
	BasisSet_t basis;
	int niters, natoms, nshells, nbasfuncs;
	int charge, nshellpairs, mat_size;
	
	// Shell quartet screening 
	double shell_screen_tol2;
	double *sp_scr_vals;  // Screening values (upper bound) of each shell pair
	
	// Integrals
	SIMINT_t simint;
	int *shell_bf_sind;   // Index of the first basis function of each shell
	int *shell_bf_num;    // Number of basis function in each shell
	
	// Matrices
	double *Hcore_mat;    // Core Hamiltonian matrix
	double *Ovlp_mat;     // Overlap matrix
	double *Fock_mat;     // Fock matrix
	double *Dens_mat;     // Density matrix
	double *J_mat;        // Coulomb matrix
	double *K_mat;        // Exchange matrix
	double *X_mat;        // Exchange matrix
	double *prev_Fock;    // Previous Fock matrices, for DIIS
	
	// Statistic 
	double mem_size;
};

typedef struct TinySCF_struct* TinySCF_t;

// Initialize TinySCF with a Cartesian basis set file (.gbs format), a molecule 
// coordinate file and the number of SCF iterations
void init_TinySCF(TinySCF_t TinySCF, char *bas_fname, char *xyz_fname, const int niters);

// Compute core Hamiltonian and overlap matrix
void TinySCF_compute_Hcore_Ovlp_mat(TinySCF_t TinySCF);

// Destroy TinySCF 
void free_TinySCF(TinySCF_t TinySCF);

#endif
