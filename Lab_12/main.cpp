#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;

double funF(double x)
{
    return log(x*x*x + 3*x*x + x + 0.1)*sin( 18*x );
}

void alloc_matrix(double ***A,int row, int col)
{
    *A = new double*[row];
    for(int i = 0; i < row; i++)
        (*A)[i]=new double[col];
}

void destroy(double **A, double row)
{
    for(int i=0;i<row;i++)
        delete [] A[i];
    
    delete [] A;
}

double method(int version, double hw,int N)
{
    double suma=0;
    double x;
    //Simpson
    if(version==1)
    {
        for(int i=0;i<(N/2)-1;i++)
            suma+=(hw / 3)*( funF(( 2*i) * hw) + 4 * funF( (2*i+1) * hw)+funF( (2*i+2) *hw));
        
    }
    //Milne
    else if(version==2)
    {
        for(int i=0;i<(N/4)-1;i++)
            suma+=((4*hw)/90)*(7*funF((4*i)*hw)+32*funF((4*i+1)*hw)+12*funF((2*i+2)*hw)+32*funF((2*i+3)*hw)+7*funF((2*i+4)*hw));
    }
    return suma;
}

void solve(const string& name, int version)
{
    ofstream out;
    out.open(name);

    double a = 0;
    double b = 1;

    double hw;
    int n=8;
    int N;

    double D[n+1][n+1];

    for(int w=0;w<=n;w++)
    {
        N=pow(2,w+version);
        hw=(b-a)/N;
        D[w][0]=method(version,hw,N);
        cout<<D[w][0]<<endl;
    }
    //ekstapolacja
    for(int k=1;k<=n;k++)
    {
        for(int w=k;w<=n;w++)
        {
            D[w][k]=(pow(4,k)*D[w][k-1]-D[w-1][k-1])/(pow(4,k)-1);
        }
    }

    cout<<endl<<endl;
    for(int i=0;i<=n;i++)
        cout<<" "<<D[i][i]<<endl;
    //Czyszczenie
    out.close();
}

int main()
{
    solve("Simpson.txt",1);
    cout<<endl<<endl;
    solve("Milne.txt",2);

    return 0;
}