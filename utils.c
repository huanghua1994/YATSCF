#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"

void copy_matrix_block(
	double *dst, const int ldd, double *src, const int lds, 
	const int nrows, const int ncols
)
{
	for (int irow = 0; irow < nrows; irow++)
		memcpy(dst + irow * ldd, src + irow * lds, DBL_SIZE * ncols);
} 

// For debug
void print_mat(double *mat, const int ldm, const int nrows, const int ncols, const char *mat_name)
{
	printf("%s:\n", mat_name);
	for (int i = 0; i < nrows; i++)
	{
		for (int j = 0; j < ncols; j++) 
		{
			int idx = i * ldm + j;
			double x = mat[idx];
			if (x >= 0.0) printf(" ");
			printf("%.7lf\t", x);
		}
		printf("\n");
	}
	printf("\n");
}
