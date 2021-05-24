#include <gsl/gsl_fft_complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double delta()
{
    return rand()/(RAND_MAX+1.0)-(0.5);
}

double funF0(double t,double T)
{
    double omega=2*M_PI/T;
    return sin(1*omega*t)+sin(2*omega*t)+sin(3*omega*t);
}

double funG(double t, double sigma)
{
    return ((1.0)/(sigma*sqrt(2*M_PI)))*exp(-((t*t)/(2*sigma*sigma)));
}

void printTAB(double *tab,int N)
{
    for(int i=0;i<N;i++)
    {
        printf("%2.3lf \n",tab[2*i]);
    }
}

void solve(int k,const char* name)
{
    FILE *plik;
    plik=fopen(name,"w");

    int N=pow(2,k);
    double T = 1.0;
    double tmax = 3*T;
    double dt = tmax/N;
    double sigma=T/20;

    double f[2*N];
    double g1[2*N];
    double g2[2*N];

    srand(time(NULL));

    for(int i=0;i<N;i++)
    {
        double t=dt*i;
        f[2*i]=funF0(t,T)+delta();
        g1[2*i]=funG(t,sigma);
        g2[2*i]=funG(t,sigma);
        fprintf(plik,"%4.5lf %4.5lf \n",t,f[2*i]);
    }
    fprintf(plik,"\n\n");
    //Transformata
    gsl_fft_complex_radix2_forward(f,N,1); 
    gsl_fft_complex_radix2_forward(g1,N,1);
    //Transformata odwrotna
    gsl_fft_complex_radix2_backward(g2,N,1);
    for(int i=0;i<N;i++)
    {
        f[2*i]=f[2*i]*(g1[2*i] + g2[2*i]);
    }
    //Transformata odwrotna na splocie
    gsl_fft_complex_radix2_backward(f,N,1);
    //Szukanie f max
    double maxi=0;
    for(int i=0;i<N;i++)
    {
        if(f[2*i]>maxi)
            maxi=f[2*i];
    }
    printf("%lf \n",maxi);
    for(int i=0;i<N;i++)
    {
        double t=dt*i;
        fprintf(plik,"%4.5lf %4.5lf \n",t,f[2*i]*2.5/maxi);
    } 
    fclose(plik); 
}

int main()
{
    solve(8, "k8.dat");
    solve(10,"k10.dat");
    solve(12,"k12.dat");
    return 0;
}