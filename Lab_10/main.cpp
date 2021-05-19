#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

double fun(double x, double y)
{
    return (5.0/2.0)*pow((x*x-y),2)+pow(1-x,2);
}

double pochodnaX(double x,double y, double dx)
{
    return (fun(x+dx,y)-fun(x-dx,y))/(2*dx);
}

double pochodnaY(double x,double y, double dy)
{
   return (fun(x,y+dy)-fun(x,y-dy))/(2*dy); 
}

double length(double x,double y)
{
    return sqrt(x*x+y*y);
}

void solve(double eps,const string& name)
{
    ofstream out;
    out.open(name);

    double x,y;
    double x0=-0.75;
    double y0=1.75;

    double h=0.1;

    for(int i=1;i<=1000;i++)
    {        
        x=x0-h*pochodnaX(x0,y0,0.0001);
        y=y0-h*pochodnaY(x0,y0,0.0001);
        out<<x<<" "<<y<<endl;
        if(length(x-x0,y-y0)<eps){
            break;
        }
        x0=x;
        y0=y;
    }
    out.close();
}

int main()
{

    solve(0.01,"eps1.dat");
    solve(0.001,"eps2.dat");
    return 0;
}