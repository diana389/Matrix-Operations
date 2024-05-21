/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>
// #include <atlas/cblas.h>

/*
 * Add your BLAS implementation here
 */
double *my_solver(int N, double *A, double *B)
{	
	// C=(At×B+B×A)×Bt

    /************ compute At×B ************/

    double *AtxB = (double*)calloc(N * N, sizeof(double));

    // Since A is upper triangular, we use cblas_dtrmm
    // Copy B to AtxB to use as the output buffer
    cblas_dcopy(N * N, B, 1, AtxB, 1);
    // Compute AtxB = At * B
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, 1.0, A, N, AtxB, N);

    /************ compute B×A ************/

    double *BxA = (double*)calloc(N * N, sizeof(double));

    // Copy B to BxA to use as the output buffer
    cblas_dcopy(N * N, B, 1, BxA, 1); 
    // Compute BxA = B * A
    cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, A, N, BxA, N);

    /************ compute sum = At×B + B×A ************/

    double *sum = (double*)calloc(N * N, sizeof(double));

    // sum = AtxB (initialize sum with the contents of AtxB)
    cblas_dcopy(N * N, AtxB, 1, sum, 1);
    // sum = sum + BxA
    cblas_daxpy(N * N, 1.0, BxA, 1, sum, 1);

    /************ compute C = sum × Bt ************/

    double *C = (double*)calloc(N * N, sizeof(double));
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, sum, N, B, N, 0.0, C, N);

    // Free memory
    free(AtxB);
    free(BxA);
    free(sum);

	return C;
}