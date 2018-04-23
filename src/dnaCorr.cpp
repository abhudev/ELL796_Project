#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <random>
#include <time.h>

typedef long long ll;

using namespace std;

// name of input file is command line argument
// window size - 50 kb
// increment - 10 kb
FILE* fin;
FILE* fout;
FILE* fout2;
FILE* fout3;
ll A, T, G, C;
ll count;
ll window;
ll windowSize;
ll increment;
ll gcount;
ll base;
ll *arr;
double *sums;
double *windowValues;
void initialiseNuc(int x1, int x2, int x3, int x4);
double cgRet();
double currentVal(ll i);

/*
	Command line arguments - 
	1 - input file
	2 - window size
	3 - increment
*/

int main(int argc, char const *argv[])
{
	fin = fopen(argv[1], "r");
	windowSize = atoll(argv[2]);
	increment = atoll(argv[3]);

	fout = fopen("GC_output.txt", "w");
	fout2 = fopen("Corr_output.txt", "w");
	fout3 = fopen("test.txt", "w");
	initialiseNuc(0, 0, 0, 0);


	window = 0;
	count = 0;

	// arr = new ll[windowSize];

	char a;

	double gcSkew;
	double corr;

	// gcount = 0;
	// base = 0;

	while(!feof(fin))
	{
		fprintf(fout3, "1\n");
		fscanf(fin, "%c", &a);
		if(feof(fin))
			break;
		switch(a)
		{
			case 'G':
					G++;
					// arr[count] = 1;
					break;
			case 'C':
					C++;	
					// arr[count] = -1;
					break;
			default:
					// arr[count] = -1;
					break;
		}
		if(a != '\n')
		{
			count++;
			gcount++;
		}
		if(count == windowSize)
		{
			gcSkew = (double)(C - G)/(double)(C + G);
			// corr = cgRet();
			// fprintf(fout2, "%lld\t%f\n", window, corr);
			fprintf(fout, "%lld\t%f\n", window, gcSkew);
			initialiseNuc(0, 0, 0, 0);
			window++;
			count = 0;
			fseek(fin, (1+increment-windowSize)*sizeof(char), SEEK_CUR);
			gcount = gcount + 1 + increment - windowSize;
			// base = gcount - 1;
		}


	}

	//------------------------GC skew ends, corr starts------------------------//
	// Try to use efficient libraries for complex numbers etc
	// Parallelize? Only if there are no loop-carried dependencies

	fseek(fin, 0, SEEK_SET);
	gcount = 0; // Total number of nucleotides
	while(!feof(fin))
	{
		fscanf(fin, "%c", &a);
		if(feof(fin))
			break;
		if(a == '\n'||a==' '||a=='\t')
			continue;
		gcount++;
	}

	arr = new ll[gcount];
	ll i, j, k, l;
	// first populate the array
	fseek(fin, 0, SEEK_SET);
	i = 0;
	while(!feof(fin))
	{
		fscanf(fin, "%c", &a);
		if(feof(fin))
			break;
		if(a == '\n'||a==' '||a=='\t')
			continue;
		if(a == 'G')
			arr[i] = 1;
		else
			arr[i] = -1;
		i++;
	}
	sums = new double[gcount];
	// the 'k' in N-k
	windowValues = new double[window];
	for(k = 1; k < windowSize-1; k++)
	{
		for(i = 0; i < gcount - k; i++)
		{
			int z = arr[i]  + arr[i+k];
			switch(z)
			{
				case 0:
						sums[i] = -1;
						break;
				case 2:
						sums[i] = 1;
						break;
			}
		}
		for(j = 0; j < window; j++)
		{
			double tmpv = 0;
			for(l = increment*j; l < increment*j + windowSize - k; l++)
			{
				if(l >= gcount - k)
					break;
				tmpv += sums[l];
			}
			if(tmpv >= 0)
				windowValues[j] += tmpv / (windowSize - k); 
			else
				windowValues[j] -= tmpv / (windowSize - k); 				
		}
	}

	for(i = 0; i < window; i++)
	{
		windowValues[i] /= (windowSize - 1);
		fprintf(fout2, "%lld\t%f\n", i+1, windowValues[i]);
	}
	printf("A = %lld, T = %lld, G = %lld, C = %lld, gcount = %lld, i = %lld, window = %lld\n", A, T, G, C, gcount, i, window);
	return 0;
}

void initialiseNuc(int x1, int x2, int x3, int x4)
{
	A = x1;
	T = x2;
	G = x3;
	C = x4;
}

double cgRet()
{
	// N - windowSize
	fprintf(fout3, "cgRet\n");
	double val = 0;
	ll i;
	// use i+1 - arrays are ZERO-INDEXED
	for(i = 1; i <= windowSize - 1; i++)
	{
		fprintf(fout3, "2\n");
		val += abs(currentVal(i));
	}
	val /= (double)(windowSize - 1);
	return val;
}

double currentVal(ll i)
{
	fprintf(fout3, "currentVal\n");
	ll j;
	double v = 0;
	fprintf(fout3, "%lld %lld %lld\n", j, i, windowSize);
	for(j = 0; j <= windowSize-i-1; j++)
	{
		fprintf(fout3, "3 %lld %lld\n", j, windowSize-i-1);
		v += (double)(arr[j]*arr[j+i]);
	}
	v /= (double)(windowSize-i);
	return v;
}