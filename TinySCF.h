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
	int niters, natoms, nshells, nbasfuncs, n_occ;
	int charge, nshellpairs, mat_size;
	
	// SCF iteration info
	int iter;
	double nuc_energy, energy;   // Nuclear energy and total energy
	double energy_delta_tol;     // SCF termination criteria for energy change
	
	// Shell quartet screening 
	double shell_scrtol2, max_scrval;
	double *sp_scrval;    // Square of screening values (upper bound) of each shell pair
	int    *uniq_sp_lid;  // Left shell id of all unique shell pairs
	int    *uniq_sp_rid;  // Right shell id of all unique shell pairs
	int    num_uniq_sp;   // Number of unique shell pairs (== nshells * (nshells+1) / 2)
	
	// ERIs
	SIMINT_t simint;
	int *shell_bf_sind;   // Index of the first basis function of each shell
	int *shell_bf_num;    // Number of basis function in each shell
	
	// Matrices and temporary arrays in SCF
	double *Hcore_mat;    // Core Hamiltonian matrix
	double *S_mat;        // Overlap matrix, no longer needed once X_mat is obtained
	double *F_mat;        // Fock matrix
	double *D_mat;        // Density matrix
	double *J_mat;        // Coulomb matrix
	double *K_mat;        // Exchange matrix
	double *X_mat;        // Basis transformation matrix
	double *Cocc_mat;     // Temporary matrix for building density matrix
	double *eigval;       // Eigenvalues for building density matrix
	int    *ev_idx;       // Index of eigenvalues, for sorting
	
	// Matrices and arrays for DIIS
	double *F0_mat;       // Previous Fock matrices
	double *R_mat;        // Intermediate matrix
	double *B_mat;        // Linear system coefficient matrix in DIIS
	double *DIIS_rhs;     // Linear system right-hand-side vector in DIIS
	int    *DIIS_ipiv;    // Permutation info for DGESV in DIIS
	
	// Statistic 
	double mem_size, init_time, S_Hcore_time, shell_scr_time;
};

typedef struct TinySCF_struct* TinySCF_t;

// Initialize TinySCF with a Cartesian basis set file (.gbs format), a molecule 
// coordinate file and the number of SCF iterations
void init_TinySCF(TinySCF_t TinySCF, char *bas_fname, char *xyz_fname, const int niters);

// Compute core Hamiltonian and overlap matrix
void TinySCF_compute_Hcore_Ovlp_mat(TinySCF_t TinySCF);

// Compute the screening values of each shell quartet and the unique shell pairs
// that survive screening using Schwarz inequality
void TinySCF_compute_sq_Schwarz_scrvals(TinySCF_t TinySCF);

// Perform SCF calculation
void TinySCF_do_SCF(TinySCF_t TinySCF);

// Destroy TinySCF 
void free_TinySCF(TinySCF_t TinySCF);

#endif
