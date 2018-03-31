#include <iostream>
#include <fstream>
// #include <bits/stdc++.h>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

typedef long long ll;

int main(int argc, char const *argv[])
{
	if(argc != 4)
	{
		printf("Usage: ./auto <input_seq> <window_size> <corr_plot>\n");
		exit(1);
	}

	// Window
	int win = atoi(argv[2]);

	// Vector sequence - G is 1, C is 
	vector<int> seq;

	char nuc;	// Input nucleotide, push in sequence

	// Read from file
	ifstream fin;
	fin.open(argv[1], ios::in);
	// Read first line of the file
	string line;
	getline(fin, line);
	cout<<"Reading sequence\n";
	// getline(fin, line); // Clear the buffer

	while(fin>>nuc)
	{
		if(nuc != 'A' && nuc != 'G' && nuc != 'C' && nuc != 'T')
			continue;
		// cout<<nuc;
		if(nuc == 'G')
			seq.push_back(1);
		else
			seq.push_back(-1);
	}

	fin.close();

	int N = seq.size();
	cout<<"Length of sequence is "<<N<<endl;

	// Compute number of windows - considering only 'whole' windows
	int num_win = N - win + 1;

	// Hold the current (win-1) correlation values, for the current window

	// Array to hold the cumulative GC skew values
	vector<float> tmp_corr(num_win);
	vector<float> corr(num_win, 0);

	cout<<"Start autocorrelation computation...\n";
	// Compute autocorrelation
	for(int k = 1; k <= win-1; ++k)
	{		
		printf("Step %d of %d\n", k, win-1);
		// IMPORTANT - division is costly, do only once
		float DC = 1.0/(win - k);	// Dividing constant
		// Zero-initialize tmp_corr
		for(int cit = 0; cit < num_win; cit++)
			tmp_corr[cit] = 0;

		for(int j = 0; j <= N-k-1; ++j)
		{
			// Unique element adding to the correlation
			float corr_element = seq[j] * seq[j+k];
			// float corr_element = seq[j] * seq[j+k] * DC;
			// TODO - fix this line
			// Fix lower limit
			int lower_i = j+k-win+1;
			if(lower_i < 0)
				lower_i = 0;
			// Fix upper limit
			int upper_i = j;
			if(upper_i > N-win)
				upper_i = N-win;
			for(int i = lower_i; i <= upper_i; ++i)
			{
				tmp_corr[i] += corr_element;
				// corr[i] += corr_element;				
			}
		}

		for(int i = 0; i < num_win; ++i)
			corr[i] += (fabs(tmp_corr[i]*DC));
		// printf("\r");
	}

	// Divide finally by win-1
	for(int i = 0; i < num_win; i++)
		corr[i] /= (win-1);
	cout<<"\nEnd autocorrelation computation\n";
	
	// Write cumulative GC to file
	FILE* fout;
	fout = fopen(argv[3], "w");	
	for(int i = 0; i < num_win; ++i)
	{		
		fprintf(fout, "%f\n", corr[i]);
	}	
	fclose(fout);
	return 0;
}