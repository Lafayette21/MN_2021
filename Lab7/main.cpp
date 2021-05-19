#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#define n 5

using std::cout;
using std::endl;

std::fstream plik;

double fun(double x)
{
    return exp(-(x*x));
}

double* rownoodlegle(double *xm,int k,double x_min,double x_max)
{
    double h=(x_max-x_min)/k;
    for(int i=0;i<n+1;i++)
    {
        xm[i]=x_min+i*h;
    }
    return xm;
}

double* Czeb_wezly(double *xm,int k,double x_min,double x_max)
{
    for(int i=0;i<k;i++)
    {
        xm[i]=0.5*((x_max-x_min)*cos(M_PI*((2*i+1)/(2*k+2)))+(x_min+x_max));
    }
    return xm;
}
double* wart_wezly(double *ym,double *xm,double k,double (*f)(double x))
{
    for(int i=0;i<k+1;i++)
    {
        ym[i]=f(xm[i]);
    }
    return ym;
}

double W_wielomian(double x,int p,double *xm,double *ym)
{
    double w=0;
    double ilo=1;
    for(int j=0;j<=p+1;j++)
    {
        ilo=1;
        for(int k=0;k<=p+1;k++)
        {
            if(k!=j) 
                ilo*=(x-xm[k])/xm[j]-xm[k];
        }
        w+=ym[j]*ilo;
    }
    return w;
}

void print(double *tab,int k)
{
    for(int i=0;i<k;i++)
        cout<<tab[i]<<" ";
    cout<<endl;
}

int main()
{
    double x_min=-5.0;
    double x_max=5.0;

    //tablica wezlow rownoodleglych
    double *xm=new double[n+1];
    xm=rownoodlegle(xm,n,x_min,x_max);
    //print(xm,n);
    

    /*Wypelnienie zgodnie z zerami wielomianow Czebyszewa
    xm=Czeb_wezly(xm,n,x_min,x_max);*/


    //Wartości funkcji w wezlach
    double *ym=new double[n+1];
    ym=wart_wezly(ym,xm,n,fun);
    //print(ym,n);
    
    //Wielomian interpolacyjny
    plik.open("zad1.dat",std::ios::out | std::ios::app);
    double x=x_min;
    do
    {
        plik<<x<<"   "<<W_wielomian(x,n,xm,ym)<<endl;
        x+=0.01; 

    }while(x<=x_max);
    plik<<endl<<endl;

    plik.close();
    delete [] xm;
    delete [] ym;
    return 0;
}