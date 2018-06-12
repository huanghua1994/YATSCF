#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>

#include <mkl.h>

#include "utils.h"
#include "TinySCF.h"
#include "build_density.h"

static void quickSort(double *eigval, int *ev_idx, int l, int r)
{
	int i = l, j = r, iswap;
	double mid = eigval[(i + j) / 2], dswap;
	while (i <= j)
	{
		while (eigval[i] < mid) i++;
		while (eigval[j] > mid) j--;
		if (i <= j)
		{
			iswap     = ev_idx[i];
			ev_idx[i] = ev_idx[j];
			ev_idx[j] = ev_idx[i];
			
			dswap     = eigval[i];
			eigval[i] = eigval[j];
			eigval[j] = dswap;
			
			i++;  j--;
		}
	}
	if (i < r) quickSort(eigval, ev_idx, i, r);
	if (j > l) quickSort(eigval, ev_idx, l, j);
}

void TinySCF_build_DenMat(TinySCF_t TinySCF)
{
	double *F_mat    = TinySCF->F_mat;
	double *D_mat    = TinySCF->D_mat;
	double *X_mat    = TinySCF->X_mat;
	double *tmp_mat  = TinySCF->tmp_mat;
	double *Cocc_mat = TinySCF->Cocc_mat;
	double *eigval   = TinySCF->eigval;
	int    *ev_idx   = TinySCF->ev_idx;
	int    nbf       = TinySCF->nbasfuncs;
	int    n_occ     = TinySCF->n_occ;
	
	// Notice: here F_mat is already = X^T * F * X
	memcpy(tmp_mat, F_mat, DBL_SIZE * TinySCF->mat_size);
	
	// Diagonalize F = C0^T * epsilon * C0, and C = X * C0 
	// [C0, E] = eig(F1), now C0 is stored in tmp_mat
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nbf, tmp_mat, nbf, eigval);  // tmp_mat will be overwritten by eigenvectors
	// C = X * C0, now C is stored in D_mat
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, X_mat, nbf, tmp_mat, nbf, 0.0, D_mat, nbf);
	
	// Form the C_occ with eigenvectors corresponding to n_occ smallest eigenvalues
	for (int i = 0; i < nbf; i++) ev_idx[i] = i;
	quickSort(eigval, ev_idx, 0, nbf - 1);
	for (int j = 0; j < n_occ; j++)
		for (int i = 0; i < nbf; i++)
			Cocc_mat[i * n_occ + j] = D_mat[i * nbf + ev_idx[j]];
	
	// D = C_occ * C_occ^T
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, n_occ, 
				1.0, Cocc_mat, n_occ, Cocc_mat, n_occ, 0.0, D_mat, nbf);
}

static double matrix_trace(double *mat, const int ldm, const int nrows)
{
	double res = 0.0;
	for (int i = 0; i < nrows; i++)
		res += mat[i * ldm + i];
	return res;
}

void TinySCF_build_DenMat_Purif(TinySCF_t TinySCF, int *purif_iter)
{
	double *F_mat    = TinySCF->F_mat;
	double *D_mat    = TinySCF->D_mat;
	double *D2_mat   = TinySCF->D2_mat;
	double *D3_mat   = TinySCF->D3_mat;
	double *X_mat    = TinySCF->X_mat;
	double *tmp_mat  = TinySCF->tmp_mat;
	int    nbf       = TinySCF->nbasfuncs;
	int    n_occ     = TinySCF->n_occ;
	
	// Gerschgorin's formula to estimate eigenvalue range
	double Hmax = -DBL_MAX, Hmin = DBL_MAX;
	for (int i = 0; i < nbf; i++)
	{
		double row_abs_sum = 0.0;
		double *row_ptr = F_mat + i * nbf;
		double Fii, Hmin0, Hmax0;
		for (int j = 0; j < nbf; j++)
			row_abs_sum += fabs(row_ptr[j]);
		Fii = row_ptr[i];
		row_abs_sum -= fabs(Fii);
		Hmin0 = Fii - row_abs_sum;
		Hmax0 = Fii + row_abs_sum;
		if (Hmin0 < Hmin) Hmin = Hmin0;
		if (Hmax0 > Hmax) Hmax = Hmax0;
	}
	
	// Generate initial guess
	double mu_bar  = matrix_trace(F_mat, nbf, nbf) / nbf;
	double lambda0 = (double) n_occ / (Hmax - mu_bar);
	double lambda1 = (double) (nbf - n_occ) / (mu_bar - Hmin);
	double lambda  = lambda0 < lambda1 ? lambda0 : lambda1;
	double coef_I  = (lambda * mu_bar + (double) n_occ) / (double) nbf;
	double coef_F  = lambda / (double) nbf;
	#pragma simd
	for (int i = 0; i < nbf * nbf; i++)
		D_mat[i] = -coef_F * F_mat[i];
	#pragma simd
	for (int i = 0; i < nbf; i++)
		D_mat[i * nbf + i] += coef_I;
	
	// Canonical Purification iterations
	int iter = 0;
	for (iter = 0; iter < MAX_PURIF_ITER; iter++)
	{
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
					1.0, D_mat, nbf, D_mat, nbf, 0.0, D2_mat, nbf);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
					1.0, D_mat, nbf, D2_mat, nbf, 0.0, D3_mat, nbf);
		
		double c, tr0, tr1;
		for (int i = 0; i < nbf * nbf; i++)
			tmp_mat[i] = D2_mat[i] - D3_mat[i];
		tr0 = matrix_trace(tmp_mat, nbf, nbf);
		for (int i = 0; i < nbf * nbf; i++)
			tmp_mat[i] = D_mat[i] - D2_mat[i];
		tr1 = matrix_trace(tmp_mat, nbf, nbf);
		c  = tr0 / tr1;
		
		if (c <= 0.5)
		{
			double c0 = 1.0 - 2.0 * c;
			double c1 = 1.0 + c;
			double c2 = 1.0 / (1.0 - c);
			#pragma simd
			for (int i = 0; i < nbf * nbf; i++)
				D_mat[i] = (c0 * D_mat[i] + c1 * D2_mat[i] - D3_mat[i]) * c2;
		} else {
			double c1 = 1.0 + c;
			double c2 = 1.0 / c;
			#pragma simd
			for (int i = 0; i < nbf * nbf; i++)
				D_mat[i] = (c1 * D2_mat[i] - D3_mat[i]) * c2;
		}
		
		// Compute the Frobenius of D - D^2 
		double err_norm = 0.0;
		for (int i = 0; i < nbf * nbf; i++)
		{
			double diff = D_mat[i] - D2_mat[i];
			err_norm += diff * diff;
		}
		err_norm = sqrt(err_norm);
		
		if (err_norm < 1e-11)   break;
		if ((c < 0) || (c > 1)) break;
	}
	
	// D = X * D * X^T
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, X_mat, nbf, D_mat, nbf, 0.0, tmp_mat, nbf);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 
				1.0, tmp_mat, nbf, X_mat, nbf, 0.0, D_mat, nbf);			
	
	*purif_iter = iter + 1;
}
