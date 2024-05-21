/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {

	// C=(At×B+B×A)×Bt

	double register *pA, *pB, *pAt, *pBt, *pAtxB, *pBxA, *pSum, *orig_pAt, *orig_pB, *orig_pSum;
	int register i, j, k;
	double register suma;

	/************ transpose matrix A ************/

	double *At = (double *)calloc(N * N, sizeof(double));
	
	for (j = 0; j < N; j++)
	{
		pAt = &At[j * N];
		pA = &A[j];

		for (i = 0; i <= j; i++)
		{
			*pAt = *pA;
			pAt++;
			pA += N;
		}
	}

	/************ compute At×B ************/

	double *AtxB = (double *)calloc(N * N, sizeof(double));

	for (i = 0; i < N; i++)
	{
		orig_pAt = &At[i * N];

		for (j = 0; j < N; j++)
		{
			pAt = orig_pAt;
			pB = &B[j];

			suma = 0.0;

			for (k = 0; k <= i; k++)
			{
				suma += *pAt * *pB;
				pAt++;
				pB += N;
			}

			AtxB[i * N + j] = suma;
		}
	}

	/************ compute B×A ************/

	double *BxA = (double *)calloc(N * N, sizeof(double));

	for (i = 0; i < N; i++)
	{
		orig_pB = &B[i * N];

		for (j = 0; j < N; j++)
		{
			pB = orig_pB;
			pA = &A[j];

			suma = 0.0;

			for (k = 0; k <= max(i, j); k++)
			{
				suma += *pB * *pA;
				pB++;
				pA += N;
			}

			BxA[i * N + j] = suma;
		}
	}

	/************ compute At×B + B×A ************/

	double *sum = (double *)calloc(N * N, sizeof(double));

	for (i = 0; i < N; i++)
	{
		pAtxB = &AtxB[i * N];
		pBxA = &BxA[i * N];

		for (j = 0; j < N; j++)
		{
			sum[i * N + j] = *pAtxB + *pBxA;
			pAtxB++;
			pBxA++;
		}
	}

	/************ transpose matrix B ************/

	double *Bt = (double *)calloc(N * N, sizeof(double));

	pB = B;

	for (i = 0; i < N; i++)
	{
		pBt = &Bt[i];

		for (j = 0; j < N; j++)
		{
			*pBt = *pB;
			pB++;
			pBt += N;
		}
	}

	/************ compute C = (At×B + B×A)×Bt ************/

	double *C = (double *)calloc(N * N, sizeof(double));

	for (i = 0; i < N; i++)
	{
		orig_pSum = &sum[i * N];

		for (j = 0; j < N; j++)
		{
			pSum = orig_pSum;
			pBt = &Bt[j];

			suma = 0.0;

			for (k = 0; k < N; k++)
			{
				suma += *pSum * *pBt;
				pSum++;
				pBt += N;
			}

			C[i * N + j] = suma;
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
