#include <iostream>
#include <fstream>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <complex>

using namespace std;

typedef long long ll;

int main(int argc, char const *argv[])
{
	if(argc != 5)
	{
		printf("Usage: ./iCorr <input_seq> <window_size> <corr_plot> <slide_dist>\n");
		exit(1);
	}

	// Window
	int win = atoi(argv[2]);
	int slide = atoi(argv[4]);

	// Vector sequence - G is 1, C is 
	vector<complex<float> > seq;

	char nuc;	// Input nucleotide, push in sequence
	// Vector of sequence bits

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
		if(nuc == 'A')
		{
			complex<float> cnum(0, 1);
			seq.push_back(cnum);
		}
		else if(nuc == 'G')
		{
			complex<float> cnum(1, 0);
			seq.push_back(cnum);
		}
		else if(nuc == 'T')
		{
			complex<float> cnum(0, -1);
			seq.push_back(cnum);
		}
		else if(nuc == 'C')
		{
			complex<float> cnum(-1, 0);
			seq.push_back(cnum);
		}
	}

	fin.close();

	int N = seq.size();
	cout<<"Length of sequence is "<<N<<endl;

	// Compute number of windows - considering only 'whole' windows
	int num_win = N - win + 1;

	// Hold the current (win-1) correlation values, for the current window

	// Array to hold the corr values - tmp and final
	vector<complex<float> > tmp_corr(num_win);
	vector<float> corr(num_win, 0);

	cout<<"Start autocorrelation computation...\n";
	// Compute autocorrelation
	for(int k = 1; k <= win-1; ++k)
	{		
		printf("Step %d of %d\n", k, win-1);
		// IMPORTANT - division is costly, do only once
		float DC = 1.0/(win - k);	// Dividing constant
		// Zero-initialize tmp_corr
		for(int cit = 0; cit < num_win; cit += slide)
			tmp_corr[cit] = 0;

		for(int j = 0; j <= N-k-1; ++j)
		{
			// Unique element adding to the correlation
			complex<float> corr_element = seq[j] * seq[j+k];
			// float corr_element = seq[j] * seq[j+k] * DC;
			// TODO - fix this line
			// Fix lower limit
			int lower_i = j+k-win+1;
			if(lower_i < 0)
				lower_i = 0;
			int rem = lower_i % slide;
			if(rem != 0)
				lower_i = lower_i - rem + slide;
			// Fix upper limit
			int upper_i = j;
			if(upper_i > N-win)
				upper_i = N-win;
			for(int i = lower_i; i <= upper_i; i += slide)
			{
				tmp_corr[i] += corr_element;
				// corr[i] += corr_element;				
			}
		}

		for(int i = 0; i < num_win; i += slide)
			corr[i] += (std::abs(tmp_corr[i]*DC));
		// printf("\r");
	}

	// Divide finally by win-1
	for(int i = 0; i < num_win; i += slide)
		corr[i] /= (win-1);
	cout<<"\nEnd autocorrelation computation\n";
	
	// Write cumulative GC to file
	FILE* fout;
	fout = fopen(argv[3], "w");	
	for(int i = 0; i < num_win; i += slide)
	{		
		fprintf(fout, "%d %f\n", i/slide, corr[i]);
	}	
	fclose(fout);
	return 0;
}