#include "blas_sparse.h"
#include <stdio.h>

main() {
	int N = 4;
	int nz = 6;
	double val[] = {1.1, 2.2, 2.4, 3.3, 4.1, 4.4};
	int indx[] = {0, 1, 1, 2, 3, 3};
	int jndx[] = {0, 1, 3, 2, 0, 3};
	double x[] = {1.0, 1.0, 1.0, 1.0};
	double y[] = {0.0, 0.0, 0.0, 0.0};

	blas_sparse_matrix A;
	int i;
	double alpha = 1.0;

	A = BLAS_duscr_begin(N, N);
	printf("Check\n");
	for (i = 0; i < nz; i++)
		BLAS_duscr_insert_entry(A, val[i], indx[i], jndx[i]);

	printf("Check\n");
	BLAS_duscr_end(A);
	printf("Check\n");

	BLAS_dusmv(blas_trans, alpha, A, x, 1, y, 1);
	printf("Check\n");

	BLAS_usds(A);
	printf("Check\n");

	for (i = 0; i < N; i++)
		printf("%lf ", y[i]);

	printf("\n");
}
