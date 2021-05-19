#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;

#define N 5
#define IT_MAX 30

double licz_r(double *a,double *b,int n, double xj)
{
    b[n]=0;
    for(int i=n-1;i>=0;i--)
    {
        b[i]=a[i+1]+xj*b[i+1];
    }
    return a[0]+xj*b[0];
}

std::ofstream out;

int main()
{
    double a[N+1]={240,-196,-92,33,14,1};
    double b[N+1];
    double c[N];
    out.open("wyniki.txt");
    //Algorytm
    for(int L=1;L<=N;L++)
    {
        int n=N-L+1;
        double x0=0;
        for(int it=1;it<IT_MAX;it++)
        {
            double Rj=licz_r(a,b,n,x0);
            double Rjp=licz_r(b,c,n-1,x0);

            double x1=x0-Rj/Rjp;
            out<<L<<" "<<it<<"   "<<x1<<"   "<<Rj<<"   "<<Rjp<<endl;

            if(fabs(x1-x0)<pow(10,-7)) 
                break;
            x0=x1;
        }
        out<<endl;
        for(int i=0;i<=n+1;i++)
            a[i]=b[i];
    }
    out.close();
    return 0;
}