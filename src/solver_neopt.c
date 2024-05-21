/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/*
 * Add your unoptimized implementation here
 */
double *my_solver(int N, double *A, double *B)
{

	// C=(At×B+B×A)×Bt

	/************ transpose matrix A ************/

	double *At = (double *)calloc(N * N, sizeof(double));

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i <= j; i++)
		{
			At[j * N + i] = A[i * N + j];
		}
	}
	
	/************ compute At×B ************/

	double *AtxB = (double *)calloc(N * N, sizeof(double));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k <= i; k++)
			{
				AtxB[i * N + j] += At[i * N + k] * B[k * N + j];
			}
		}
	}

	/************ compute B×A ************/
	
	double *BxA = (double *)calloc(N * N, sizeof(double));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k <= max(i, j); k++)
			{
				BxA[i * N + j] += B[i * N + k] * A[k * N + j];
			}
		}
	}
	
	/************ compute sum = At×B + B×A ************/

	double *sum = (double *)calloc(N * N, sizeof(double));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			sum[i * N + j] = AtxB[i * N + j] + BxA[i * N + j];
		}
	}

	/************ transpose matrix B ************/
	
	double *Bt = (double *)calloc(N * N, sizeof(double));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Bt[j * N + i] = B[i * N + j];
		}
	}

	/************ compute C = (At×B + B×A)×Bt ************/

	double *C = (double *)calloc(N * N, sizeof(double));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				C[i * N + j] += sum[i * N + k] * Bt[k * N + j];
			}
		}
	}

	// free memory
	free(At);
	free(AtxB);
	free(BxA);
	free(sum);
	free(Bt);

	return C;
}
