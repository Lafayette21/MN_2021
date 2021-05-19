#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>

#include "nrutil.h"
#include "nrutil.c" 
#include "gaussj.c" 

using namespace std;

#define N 100 // rozmiar macierzy M: NxN

int main()
{
	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);
	float h=0.1;

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
	}
	M[2][1]=-1.0;
    b[1][1]=1.0;
	
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			cout<<M[i][j]<<" ";
		}
		cout<<endl;
	}

	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
    cout<<"Wektor b: "<<endl;
	for (int i = 1; i <= N; ++i)
        cout<<b[i][1]<<endl;
	ofstream out;
	out.open("wyniki.dat");
	for(int i=1;i<=N;i++)
		out<<b[i][1]<<endl;

	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	out.close();
	return 0;
}