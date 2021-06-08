#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void gaulag(float x[], float w[], int n, float alf)
{
    const double EPS = 3.0e-14;
    const double MAXIT = 10;
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i,its,j;
	float ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) nrerror("too many iterations in gaulag");
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((float)n))/(pp*n*p2);
	}
}
void gauher(float x[], float w[], int n)
{
    const double EPS = 3.0e-14;
    const double PIM4 = 0.7511255444649425;
    const double MAXIT = 10;
	void nrerror(char error_text[]);
	int i,its,j,m;
	double p1,p2,p3,pp,z,z1;

	m=(n+1)/2;
	for (i=1;i<=m;i++) {
		if (i == 1) {
			z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
		} else if (i == 2) {
			z -= 1.14*pow((double)n,0.426)/z;
		} else if (i == 3) {
			z=1.86*z-0.86*x[1];
		} else if (i == 4) {
			z=1.91*z-0.91*x[2];
		} else {
			z=2.0*z-x[i-2];
		}
		for (its=1;its<=MAXIT;its++) {
			p1=PIM4;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
			}
			pp=sqrt((double)2*n)*p2;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) nrerror("too many iterations in gauher");
		x[i]=z;
		x[n+1-i] = -z;
		w[i]=2.0/(pp*pp);
		w[n+1-i]=w[i];
	}
}

void gauleg(float x1, float x2, float x[], float w[], int n)
{
    const double EPS = 3.0e-11;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}

float fun1(float x)
{
    return 1.0/(x*sqrt(x*x-1));
}

float fun2her(float x)
{
    return 0.5*log(fabs(x));
}
float fun2leg(float x)
{
    return log(x)*exp(-x*x);
}

float fun3(float x)
{
    return sin(2*x)*exp(-2*x); //wagowa bezposrednio wydzielona
}

int main()
{
    //Zad1 fun1 gauleg
    float *x;
    float *w; 
    float suma,blad;
    float c1=M_PI/3;

    FILE *plik;
    plik=fopen("out.dat","w");

    for(int n=2;n<=100;n++)
    {
        suma=0;
        x=vector(1,n);
        w=vector(1,n);
        gauleg(1.0,2.0,x,w,n);
        //liczenie calki
        for(int i=1;i<n;i++)
        {
            printf("%lf %lf \n\n",x[i],w[i]);
            suma+=w[i]*fun1(x[i]);
        }
        printf(" %d nowa \n",n);
        blad=fabs(c1-suma);
        fprintf(plik,"%d %f \n",n,blad);
    }
    fprintf(plik,"\n\n");
    //Zad 2 fun2 gauher vs gauleg
    float c2=-0.8700577;
    //Gauher
    for(int n=2;n<=100;n+=2)
    {
        suma=0;
        x=vector(1,n);
        w=vector(1,n);
        gauher(x,w,n);
        //liczenie calki
        for(int i=1;i<n;i++)
        {
            suma+=w[i]*fun2her(x[i]);
        }
        blad=fabs(c2-suma);
        fprintf(plik,"%d %f \n",n,blad);
    }
    fprintf(plik,"\n\n");
    //Gauleg
    for(int n=2;n<=100;n++)
    {
        suma=0;
         x=vector(1,n);
        w=vector(1,n);
        gauleg(0.0,5.0,x,w,n);
        //liczenie calki
        for(int i=1;i<n;i++)
        {
            suma+=w[i]*fun2leg(x[i]);
        }
        blad=fabs(c2-suma);
        fprintf(plik,"%d %f \n",n,blad);
    }
    fprintf(plik,"\n\n");
    //Zad3 fun3 gaulag
    double c3=2/13;
    for(int n=2;n<=20;n++)
    {
        suma=0;
        x=vector(1,n);
        w=vector(1,n);
        gaulag(x,w,n,0.0);
        //liczenie calki
        for(int i=1;i<n;i++)
        {
            suma+=w[i]*fun3(x[i]);
        }
        blad=fabs(c3-suma);
        fprintf(plik,"%d %f \n",n,blad);
    }
    fclose(plik);
    return 0;
}