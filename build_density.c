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

static void Gerschgorin_MaxMin(double *M, int ldm, int nrows, double *max_ev, double *min_ev)
{
	double Hmax = -DBL_MAX, Hmin = DBL_MAX;
	for (int i = 0; i < nrows; i++)
	{
		double row_abs_sum = 0.0;
		double *row_ptr = M + i * ldm;
		double Fii, Hmin0, Hmax0;
		for (int j = 0; j < nrows; j++)
			row_abs_sum += fabs(row_ptr[j]);
		Fii = row_ptr[i];
		row_abs_sum -= fabs(Fii);
		Hmin0 = Fii - row_abs_sum;
		Hmax0 = Fii + row_abs_sum;
		if (Hmin0 < Hmin) Hmin = Hmin0;
		if (Hmax0 > Hmax) Hmax = Hmax0;
	}
	*max_ev = Hmax;
	*min_ev = Hmin;
}

void TinySCF_build_DenMat_Purif(TinySCF_t TinySCF, int *purif_iter)
{
	double *F_mat   = TinySCF->F_mat;
	double *D_mat   = TinySCF->D_mat;
	double *D2_mat  = TinySCF->D2_mat;
	double *D3_mat  = TinySCF->D3_mat;
	double *X_mat   = TinySCF->X_mat;
	double *tmp_mat = TinySCF->tmp_mat;
	int    nbf      = TinySCF->nbasfuncs;
	int    n_occ    = TinySCF->n_occ;
	
	// Gerschgorin's formula to estimate eigenvalue range
	double Hmax, Hmin;
	Gerschgorin_MaxMin(F_mat, nbf, nbf, &Hmax, &Hmin);
	
	// Generate initial guess
	double mu_bar  = 0.0;
	for (int i = 0; i < nbf; i++)
		mu_bar += F_mat[i * nbf + i];
	mu_bar /= (double) nbf;
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
		
		double c, tr0 = 0.0, tr1 = 0.0;
		for (int i = 0; i < nbf; i++)
		{
			tr0 += D2_mat[i * nbf + i] - D3_mat[i * nbf + i];
			tr1 += D_mat[i * nbf + i]  - D2_mat[i * nbf + i];
		}
		c = tr0 / tr1;
		
		double c0, c1, c2;
		if (c <= 0.5)
		{
			c0 = 1.0 - 2.0 * c;
			c1 = 1.0 + c;
			c2 = 1.0 / (1.0 - c);
		} else {
			c1 = 1.0 + c;
			c2 = 1.0 / c;
		}
		
		double err_norm = 0.0;
		
		#pragma omp parallel
		{
			if (c <= 0.5)
			{
				#pragma omp for reduction(+:err_norm)
				for (int i = 0; i < nbf * nbf; i++)
				{
					D_mat[i] = (c0 * D_mat[i] + c1 * D2_mat[i] - D3_mat[i]) * c2;
					double diff = D_mat[i] - D2_mat[i];
					err_norm += diff * diff;
				}
			} else {
				#pragma omp for reduction(+:err_norm)
				for (int i = 0; i < nbf * nbf; i++)
				{
					D_mat[i] = (c1 * D2_mat[i] - D3_mat[i]) * c2;
					double diff = D_mat[i] - D2_mat[i];
					err_norm += diff * diff;
				}
			}
		}
		err_norm = sqrt(err_norm);
		
		if (err_norm < PURIF_TOL)   break;
		if ((c < 0) || (c > 1)) break;
	}
	
	// D = X * D * X^T
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, X_mat, nbf, D_mat, nbf, 0.0, tmp_mat, nbf);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 
				1.0, tmp_mat, nbf, X_mat, nbf, 0.0, D_mat, nbf);			
	
	*purif_iter = iter + 1;
}

void TinySCF_build_DenMat_SP2(TinySCF_t TinySCF, int *SP2_iter)
{
	double *F_mat   = TinySCF->F_mat;
	double *D_mat   = TinySCF->D_mat;
	double *X_mat   = TinySCF->X_mat;
	double *tmp_mat = TinySCF->tmp_mat;
	int    nbf      = TinySCF->nbasfuncs;
	double n_elec   = 2.0 * TinySCF->n_occ;

	// Gerschgorin's formula to estimate eigenvalue range
	double Hmax, Hmin;
	Gerschgorin_MaxMin(F_mat, nbf, nbf, &Hmax, &Hmin);

	// Generate initial guess
	double trace_D = 0.0, inv_ev_diff = 1.0 / (Hmax - Hmin);
	for (int i = 0; i < nbf * nbf; i++)
		D_mat[i] = -F_mat[i] * inv_ev_diff;
	for (int i = 0; i < nbf; i++)
	{
		D_mat[i * nbf + i] += Hmax * inv_ev_diff;
		trace_D += D_mat[i * nbf + i];
	}

	// SP2 iterations
	int iter = 0;
	for (iter = 0; iter < MAX_SP2_ITER; iter++)
	{
		memcpy(tmp_mat, D_mat, DBL_SIZE * nbf * nbf);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
					-1.0, D_mat, nbf, D_mat, nbf, 1.0, tmp_mat, nbf);

		double trace_tmpD = 0;
		for (int i = 0; i < nbf; i++)
			trace_tmpD += tmp_mat[i * nbf + i];

		double t1 = fabs(2.0 * (trace_D - trace_tmpD) - n_elec);
		double t2 = fabs(2.0 * (trace_D + trace_tmpD) - n_elec);
		if (t1 > t2)
		{
			for (int i = 0; i < nbf * nbf; i++)
				D_mat[i] += tmp_mat[i];
			trace_D += trace_tmpD;
		} else {
			for (int i = 0; i < nbf * nbf; i++)
				D_mat[i] -= tmp_mat[i];
			trace_D -= trace_tmpD;
		}

		double idem_err = fabs(trace_tmpD);
		if (idem_err < SP2_TOL) break;
	}

	// D = X * D * X^T
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, X_mat, nbf, D_mat, nbf, 0.0, tmp_mat, nbf);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 
				1.0, tmp_mat, nbf, X_mat, nbf, 0.0, D_mat, nbf);			
	
	*SP2_iter = iter + 1;
}

void TinySCF_build_DenMat_SSNS(TinySCF_t TinySCF, int *SSNS_iter, int use_scale)
{
	double *F_mat   = TinySCF->F_mat;
	double *D_mat   = TinySCF->D_mat;
	double *D2_mat  = TinySCF->D2_mat;
	double *D3_mat  = TinySCF->D3_mat;
	double *X_mat   = TinySCF->X_mat;
	double *tmp_mat = TinySCF->tmp_mat;
	double *eigval  = TinySCF->eigval;
	int    *ev_idx  = TinySCF->ev_idx;
	int    nbf      = TinySCF->nbasfuncs;
	int    n_occ    = TinySCF->n_occ;
	
	double mu, absHmax = 0, absHmin = 9e99;
	double x, alpha_hat = 1.69770248525577;
	
	// Solve the eigen system to get mu (chemical potential) and the max/min
	// of negative and non-negative eigenvalues. In practice these values should
	// be obtained by other faster methods
	memcpy(tmp_mat, F_mat, DBL_SIZE * TinySCF->mat_size);
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nbf, tmp_mat, nbf, eigval); 
	quickSort(eigval, ev_idx, 0, nbf - 1);   // Sort the eigenvalues to ascending; we don't need ev_idx here
	mu = (eigval[n_occ - 1] + eigval[n_occ]) * 0.5;
	for (int i = 0; i < nbf; i++)   // Shift the eigenvalues for mu * I - F
	{
		eigval[i] = mu - eigval[i];
		double abs_ev = fabs(eigval[i]);
		if (abs_ev > absHmax) absHmax = abs_ev;
		if (abs_ev < absHmin) absHmin = abs_ev;
	}
	x = absHmin / absHmax;
	
	double st = get_wtime_sec();
	
	// Generate initial guess
	double inv_absHmax = 1.0 / absHmax;
	for (int i = 0; i < nbf * nbf; i++)
		D_mat[i] = -F_mat[i] * inv_absHmax;
	for (int i = 0; i < nbf; i++)
		D_mat[i * nbf + i] += mu * inv_absHmax;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, D_mat, nbf, D_mat, nbf, 0.0, D2_mat, nbf);
	
	// Purification iterations
	int iter = 0;
	for (iter = 0; iter < MAX_SP2_ITER; iter++)
	{
		double alpha = sqrt(3.0 / (1.0 + x + x * x));
		if (alpha > alpha_hat) alpha = alpha_hat;
		if (use_scale == 0) alpha = 1.0;
		
		double semi_alpha = 0.5 * alpha;
		double neg_alpha2 = -alpha * alpha;
		for (int i = 0; i < nbf * nbf; i++)
		{
			tmp_mat[i] = D_mat[i] * semi_alpha;
			D2_mat[i] *= neg_alpha2;
		}
		for (int i = 0; i < nbf; i++)
			D2_mat[i * nbf + i] += 3.0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
					1.0, tmp_mat, nbf, D2_mat, nbf, 0.0, D_mat, nbf);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
					1.0, D_mat, nbf, D_mat, nbf, 0.0, D2_mat, nbf);
		x = semi_alpha * x * (3.0 + neg_alpha2 * x * x);
		
		// realD = (D + I)/2, realD2 = D^2/4 + D/2 + I/4
		// realD2 - realD = D^2/4 - I/4, compute its Frobenius norm 
		double err_norm = 0.0;
		for (int i = 0; i < nbf; i++)
			D2_mat[i * nbf + i] -= 1.0;
		for (int i = 0; i < nbf * nbf; i++)
			err_norm += D2_mat[i] * D2_mat[i];
		err_norm /= 16.0;
		for (int i = 0; i < nbf; i++)
			D2_mat[i * nbf + i] += 1.0;
		
		if (err_norm < 1e-11) break;
	}
	
	// D = (D + I) / 2
	for (int i = 0; i < nbf * nbf; i++)
		D_mat[i] *= 0.5;
	for (int i = 0; i < nbf; i++)
		D_mat[i * nbf + i] += 0.5;
	
	// D = X * D * X^T
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 
				1.0, X_mat, nbf, D_mat, nbf, 0.0, tmp_mat, nbf);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 
				1.0, tmp_mat, nbf, X_mat, nbf, 0.0, D_mat, nbf);			
	
	double et = get_wtime_sec();
	printf("  SSNS actual time = %lf (s)\n", et - st);
	
	*SSNS_iter = iter + 1;
}

