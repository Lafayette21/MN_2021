#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c" 
#include "gaussj.c" 

#define N 7 // rozmiar macierzy M: NxN

int main()
{
	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);
	float h=0.0;

	// 	Wypelnienie macierzy M i wektora b
	for (int i = 1; i <= N; ++i)
	{
		b[i][1] = 0.0;
		for (int j = 1; j <= N; ++j)
		{
			if(i==j)
				M[i][j]=1.0;	
			else if(j==i-2)
				M[i][j]=1.0;
			else if(j==i-1)
				M[i][j]=1.0*h*h-2.0;
			else
				M[i][j] = 0.0;
		}
		h+=0.1;
	}
	M[2][1]=-1;
	b[1][1]=1;
	
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			printf("%1.2f ",M[i][j]);
		}
		printf("\n");
	}

	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	for (int i = 1; i <= N; ++i)
		printf("%g\n", b[i][1]);

	FILE *plik;	

	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}