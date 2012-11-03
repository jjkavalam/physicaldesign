/*

    Compile instructions:
        gcc myblas-test.c -I . -L . -lmyblas
        
    This test checks the implementation
    * Compressed Row Storage of a sparse matrix
    * sparse-matrix vector multiplication

*/


#include <stdio.h>

#include "myblas.h"

main() {
	int N = 4;

	double x[] = {1.0, 1.0, 1.0, 1.0};
	double y[] = {0.0, 0.0, 0.0, 0.0};

    int cols[] = {0}; 
    double vals[] = {0.0};
    
	int i;
	double alpha = 1.0;

	BLAS_matrix A = BLAS_duscr_begin(N, N);
    //BLAS_duscr_insert_row( A, i, 1, vals, cols )
    //BLAS_duscr_end(A);
	//BLAS_dusmv(blas_trans, alpha, A, x, 1, y, 1);
	BLAS_usds(A);
}
