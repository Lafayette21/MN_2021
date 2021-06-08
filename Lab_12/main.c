#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define n 8

const double a = 0.0;
const double b = 1.0;

double fun(double x)
{
	return log(pow(x, 3.0) + 3 * pow(x, 2.0) + x + 0.1) * sin(18 * x);
}

void Richardson(double** D)
{
    int i, k;
	for (i = 1; i <= n; i++)
	{
		for (k = 1; k < i + 1; k++)
        {
			D[i][k] = (pow(4, k) * D[i][k-1] - D[i-1][k-1]) / (pow(4, k) - 1);
        }
	}
}

void Simpson(double** D)
{
    int w, i;
    for (w = 0; w <= n; w++)
	{
		double h = (b - a) / pow(2.0, (w + 1));
		int N = pow(2, w + 1);

		double S = 0.0;

		for (i = 0; i <= (N/2) - 1; i++)
		{
			S += (h / 3) * (fun(a + (2 * i) * h) + 4 * fun(a + (2 * i + 1) * h)+ fun(a + (2 * i + 2) * h));
		}
		D[w][0] = S;
	}
	Richardson(D);
}

void Milne(double** D)
{
    int w, i;
    for (w = 0; w <= n; w++)
	{
		double h = (b - a) / pow(2.0, (w + 2));
		int N = pow(2, w + 2);

		double S = 0.0;

		for (i = 0; i <= (N / 4) - 1; i++)
		{
			S += ((4 * h) / 90) * (7 * fun(a + (4 * i) * h ) + 32 * fun(a + (4 * i + 1) * h) + 12 * fun(a + (4 * i + 2) * h) + 32 * fun(a + (4 * i + 3) * h) + 7 * fun(a + (4 * i + 4) * h));
		}

		D[w][0] = S;
	}

	Richardson(D);
}

void to_file(FILE* file, double **D)
{
    int i, j;
	for (i = 0; i <= n; i++)
	{
		for (j = 0; j < i+1; j++)
        {
			fprintf(file, "%.10g\t ", D[i][j]);
        }
		fprintf(file, "\n");
	}
}

int main()
{
    FILE* fp1 = fopen("Simpson.dat", "w");
    FILE* fp2 = fopen("Milne.dat", "w");

	double** D = (double**)malloc((n+1)* sizeof(double*));

	for (int i = 0; i <= n; i++)
    {
		D[i] = (double*)malloc((i+1) * sizeof(double));
    }

    Simpson(D);
	to_file(fp1, D);

    Milne(D);
	to_file(fp2, D);

	for(int i = 0; i <= n; i++)
    {
		free(D[i]);
    }
    free(D);

    fclose(fp1);
    fclose(fp2);

	return 0;
}
