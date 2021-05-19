#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>


double delta(double x,double alpha)
{
    return alpha*(rand()/(RAND_MAX+1.0)-0.5);
}


double funG(double x,double x0,double sigma,double alpha)
{
    double mian=2*sigma*sigma;
    double licz=-(pow((x-x0),2));
    return exp(licz/mian)*(1+delta(x,alpha));
}

double aproxfuG(gsl_vector *b,double x,int m)
{
    return exp(gsl_vector_get(b,0)+gsl_vector_get(b,1)*x+gsl_vector_get(b,2)*x*x+gsl_vector_get(b,3)*x*x*x);
}


void printV(gsl_vector* V,int n)
{
    for(int i=0;i<n;i++)
    {
        printf("%6.5lf ",gsl_vector_get(V,i));
    }
    printf("\n \n");
}
void printM(gsl_matrix *A,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            printf("%14.4lf ",gsl_matrix_get(A,i,j));
        printf("\n");
    }
}   

void solve(double x0,double sigma,double alpha,int N)
{
    FILE *out,*pkt;
    
    int m=4;

    gsl_matrix *G=gsl_matrix_calloc(m,m);
    gsl_vector *r=gsl_vector_calloc(m);

    gsl_vector *xw=gsl_vector_calloc(N);
    gsl_vector *gw=gsl_vector_calloc(N);
    gsl_vector *fw=gsl_vector_calloc(N);

    double xmin=-3.0*sigma+x0;
    double xmax=3.0*sigma+x0;
    double deltaX=(xmax-xmin)/(N-1);
    
    out=fopen("G.dat","a");
    pkt=fopen("pkt.dat","a");


    double x;
    for(int i=0;i<N;i++)
    {
        x=xmin+i*deltaX;
        gsl_vector_set(xw,i,x);
        gsl_vector_set(gw,i,funG(x,x0,sigma,alpha));
        gsl_vector_set(fw,i,log(gsl_vector_get(gw,i)));
        //fprintf(pkt,"%lf %lf \n",gsl_vector_get(xw,i),gsl_vector_get(gw,i));
    }
    //fprintf(pkt,"\n \n");

    for(int k=0;k<=m-1;k++)
    {
        double sumR=0;
        for(int j=0;j<=N-1;j++)
        {
            sumR+=gsl_vector_get(fw,j)*pow(gsl_vector_get(xw,j),k);
            gsl_vector_set(r,k,sumR);
        }
        for(int i=0;i<=m-1;i++)
        {
            double sumGw=0;
            for(int j=0;j<=N-1;j++)
            {
                sumGw+=pow(gsl_vector_get(xw,j),i+k);
                gsl_matrix_set(G,i,k,sumGw);
            }
        }
    }
    //Rozwiazanie HouseHolder
    gsl_linalg_HH_svx(G,r);

    printV(r,m);

    x=xmin;
    while(x<=xmax)
    {
        //fprintf(out,"%10.4lf %10.4lf \n",x,aproxfuG(r,x,m));
        x+=0.01;
    }
    //fprintf(out,"\n \n");

    fclose(out);
    fclose(pkt); 
}




int main()
{
    double x0=2.0;
    double sigma=4.0;
    //N==11 alpha==0
    solve(x0,sigma,0,11);
    //N==11 alpha ==0.5
    solve(x0,sigma,0.5,11);
    //N==101alpha == 0.5
    solve(x0,sigma,0.5,101);

    return 0;
}