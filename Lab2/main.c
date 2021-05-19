
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

#define N 4

void print(gsl_matrix *A,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%10.5lf",gsl_matrix_get(A,i,j));
        }
	printf("\n");
    }
    printf("\n\n\n");
}

double norm(gsl_matrix *A,int n)
{
    int i;
    gsl_vector *v=gsl_vector_calloc(N);
    gsl_vector *wynik=gsl_vector_calloc(N);
    for(int i=0;i<N;i++)
    {
        double suma=0.0;
        gsl_matrix_get_col(v,A,i);
        for(int j=0;j<N;j++)
        {
            suma+=fabs(gsl_vector_get(v,j));
        }
        gsl_vector_set(wynik,i,suma);
    }
    double w=gsl_vector_max(wynik);
    gsl_vector_free(v);
    gsl_vector_free(wynik);
    return w;
}

int main()
{
    //Alokacja i wypelnienie danych
    printf("Macierz A \n\n");
    gsl_permutation *p=gsl_permutation_calloc(N);
    gsl_matrix *A=gsl_matrix_calloc(N,N);
    float pomoc;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            pomoc=1.0/(i+j+2.0);
            gsl_matrix_set(A,i,j,pomoc);
        }
    }
    for(int i=0;i<N;i++)
    {
	for(int j=0;j<N;j++)
	{
	    printf("%10.6lf ",gsl_matrix_get(A,i,j));
	}
	printf("\n");
    } 
     //Kopia tempA
    gsl_matrix *tempA=gsl_matrix_calloc(N,N);
    for(int i=0;i<N;i++)
    {
    	for(int j=0;j<N;j++)
   	{
        	gsl_matrix_set(tempA,i,j,gsl_matrix_get(A,i,j));
    	}
    }
    //zad 1 Dekompozycja LU
    printf("\n \n Dekompozycja LU \n \n");
    int signum;
    gsl_linalg_LU_decomp(A,p,&signum);
    for(int i=0;i<N;i++)
    {
    	for(int j=0;j<N;j++) {
        	printf("%16.9lf",gsl_matrix_get(A,i,j));
	}
	printf("\n");
    }
    //zad 2 wyznacznik
    printf("\n Wyznacznik \n");
    double p_det=1.0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j) 
		p_det*=gsl_matrix_get(A,i,j);
        }
    }
    double det=signum*p_det;
    printf("\n %.20lf",det);
    //Zad 3 Macierz odwrotna
    gsl_matrix *Odw=gsl_matrix_calloc(N,N);
    printf("\n \n Macierz odwrotna \n \n");
    for(int i=0;i<N;i++)
    {
    	gsl_vector *b=gsl_vector_calloc(N);
	gsl_vector_set(b,i,1.0);
    	gsl_vector *x=gsl_vector_calloc(N);
    	gsl_linalg_LU_solve(A,p,b,x);
    	for(int j=0;j<N;j++)
    	{
		gsl_matrix_set(Odw,j,i,gsl_vector_get(x,j));
    	}
	gsl_vector_free(b);
	gsl_vector_free(x);
    }
    for(int i=0;i<N;i++)
    {
    	for(int j=0;j<N;j++) {
        	printf("%10.2lf",gsl_matrix_get(Odw,i,j));
        }
	printf("\n");
    }
    //Zad 4 Mnozenie macierzy odwrotnej i pierwotnej
    printf("\n \nMnozenie macierzy i odwrotnej \n");
    gsl_matrix *temp=gsl_matrix_calloc(N,N);
    for(int i=0;i<N;i++)
    {
   	 for(int j=0;j<N;j++)
   	 {
   	    	 gsl_matrix_set(temp,i,j,gsl_matrix_get(A,i,j));
    	 }
    }
    print(temp,N);
    gsl_matrix_mul_elements(temp,Odw);
    for(int i=0;i<N;i++)
    {
    	for(int j=0;j<N;j++) {
        	printf("%18.5lf",gsl_matrix_get(temp,i,j));
	}
	printf("\n");
    }	 
    gsl_matrix_free(temp);
    //Zad 5 Wskaznik uwarunkowania macierzy 
    printf("\n\n Wskaznik uwarunkowania macierzy \n\n");
    double nOdw=gsl_matrix_max(Odw);
    printf("\n %.3lf",nOdw);
    double nA=gsl_matrix_max(tempA);
    printf("\n\n %lf",nA);
    //Wynik
    printf("\n \n%10.6lf",nOdw*nA);
    

     //Zwalnianie pamieci
     gsl_matrix_free(A);
     gsl_permutation_free(p);
     gsl_matrix_free(tempA);
     gsl_matrix_free(Odw);

    return 0;
}
