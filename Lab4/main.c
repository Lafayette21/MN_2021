#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_eigen.h>

#define L 10
#define n 200
#define N 1

double dx=L/(n+1.0);

double delta(int i, int j)
{
    if(i==j) return 1.0;
    else return 0.0;
}
double ro(double x,double alfa)
{
    return 1+4*alfa*x*x;
}
double wx(int i)
{
    return -(L/2.0)+dx*(i+1);
}

void print_M(gsl_matrix *A,int k)
{
    for(int i=0;i<k;i++)
    {
        for(int j=0;j<k;j++)
        {
            printf("%.2lf ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
}

int main()
{
    FILE* plik1,*plik2;
    plik1=fopen("eval.dat","w");
    plik2=fopen("evec.dat","w");

    //Macierz A i B
    double dalfa=0.0;
    
    while(dalfa<=100)
    {
        gsl_matrix *A=gsl_matrix_calloc(n,n);
        gsl_matrix *B=gsl_matrix_calloc(n,n);
        gsl_vector *u=gsl_vector_calloc(n);
        gsl_matrix *evec=gsl_matrix_calloc(n,n);
        gsl_eigen_gensymmv_workspace *w=gsl_eigen_gensymmv_alloc(n);
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                gsl_matrix_set(A,i,j,((-delta(i,j+1)+2*delta(i,j)-delta(i,j-1))/(dx*dx)));
                gsl_matrix_set(B,i,j,(ro(wx(i),dalfa)/N)*delta(i,j));
            }
        }
        
        gsl_eigen_gensymmv(A,B,u,evec,w);
        //Sortowanie wektorow i wartosci wlasnych
        gsl_eigen_gensymmv_sort(u,evec,GSL_EIGEN_SORT_VAL_ASC);
        //plik eval.dat
        fprintf(plik1,"%.3lf ",dalfa);
        for(int i=0;i<6;i++)
            fprintf(plik1,"%7.3lf",gsl_vector_get(u,i));
        fprintf(plik1,"\n");
        //plik evec.dat
        if(dalfa==0 || dalfa==100)
        {
            for(int i=0;i<n;i++)
            {
                fprintf(plik2,"%6.2lf ", (-L/2)+dx*(i+1));
                for(int j=0;j<6;j++)
                {
                    
                    fprintf(plik2,"%6.2lf ",gsl_matrix_get(evec,i,j));
                }
                fprintf(plik2,"\n");
            }
            fprintf(plik2,"\n\n");
            
        }

        dalfa+=2.0;  
        gsl_matrix_free(A);
        gsl_matrix_free(B);
        gsl_vector_free(u);
        gsl_eigen_gensymmv_free(w);
        gsl_matrix_free(evec);
    }
    fclose(plik1);
    fclose(plik2);

    
    return 0;
}