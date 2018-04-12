#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include <mkl.h>

#include "utils.h"
#include "TinySCF.h"
#include "build_Fock.h"

void TinySCF_DIIS(TinySCF_t TinySCF)
{
	double *F_mat     = TinySCF->F_mat;
	double *F0_mat    = TinySCF->F0_mat;
	double *R_mat     = TinySCF->R_mat;
	double *B_mat     = TinySCF->B_mat;
	double *DIIS_rhs  = TinySCF->DIIS_rhs;
	int    *ipiv      = TinySCF->DIIS_ipiv;
	int mat_size      = TinySCF->mat_size;
	int mat_mem_size  = DBL_SIZE * mat_size;
	int ldB           = MAX_DIIS + 1;
	
	// F_mat is treated as a column vector in my MATLAB code
	// For performance, we treat it as a row vector here, and
	// perform necessary transpose for the R_mat
	
	// Drop the oldest F_mat and its corresponding residual vector
	for (int i = 0; i < MAX_DIIS - 1; i++)
	{
		memcpy(R_mat  + i * mat_size, R_mat  + (i + 1) * mat_size, mat_mem_size);
		memcpy(F0_mat + i * mat_size, F0_mat + (i + 1) * mat_size, mat_mem_size);
	}
	
	// Copy new F_mat and compute its residual vector
	int offset0 = (MAX_DIIS - 2) * mat_size;
	int offset1 = (MAX_DIIS - 1) * mat_size;
	memcpy(F0_mat + offset1, F_mat, mat_mem_size);
	for (int i = 0; i < mat_size; i++)
		R_mat[offset1 + i] = F0_mat[offset1 + i] - F0_mat[offset0 + i];
	
	// Start DIIS
	int ndiis = TinySCF->iter > MAX_DIIS ? MAX_DIIS : TinySCF->iter;
	if (ndiis <= 2) return;
	
	// Reconstruct B_mat and DIIS_rhs, since they will be overwritten later
	memset(B_mat, 0, DBL_SIZE * ldB * ldB);
	memset(DIIS_rhs, 0, DBL_SIZE * ldB);
	for (int i = 0; i < MAX_DIIS; i++)
	{
		B_mat[i * ldB + MAX_DIIS] = -1.0;
		B_mat[MAX_DIIS * ldB + i] =  1.0;
	}
	DIIS_rhs[ldB - 1] = 1;

	// B(1 : max_diis, 1 : max_diis) = R^T * R;
	// Here our R_mat equals to R^T since R_mat is stored in row-major style
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, MAX_DIIS, MAX_DIIS, mat_size,
				1.0, R_mat, mat_size, R_mat, mat_size, 0.0, B_mat, ldB);
	
	// Only the ndiis * ndiis block in the bottom-right corner of R^T * R
	// is non-zero, so we use this part to solve for DIIS coefficients
	int diis_spos   = MAX_DIIS - ndiis;
	int num_equ     = ndiis + 1;
	double *B_ptr   = B_mat + diis_spos * ldB + diis_spos;
	double *rhs_ptr = DIIS_rhs  + diis_spos;
	// Solve B * c = rhs with B_ptr, B_ptr will be overwritten as P * L * U, 
	// rhs_ptr will be overwritten as solution, ipiv is the permuation info
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, num_equ, 1, B_ptr, ldB, ipiv, rhs_ptr, 1);
	
	// Reconstruct F_mat, F_new = \sum_{i} c(i) * F0(:, i);
	memset(F_mat, 0, mat_mem_size);
	for (int i = diis_spos; i < MAX_DIIS; i++)
		cblas_daxpy(mat_size, DIIS_rhs[i], F0_mat + i * mat_size, 1, F_mat, 1);
}
