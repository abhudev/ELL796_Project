#include <iostream>
#include <fstream>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include "transform_seq.h"
#include <random>

using namespace std;

typedef long long ll;

int main(int argc, char const *argv[])
{
	if(argc != 7)
	{
		printf("Usage: ./ext <input_seq> <window_size> <corr_plot> <slide_dist> <kmer> <overlap>\n");
		exit(1);
	}

	// Window
	int win = atoi(argv[2]);
	int slide = atoi(argv[4]);
	int kmer = atoi(argv[5]);
	int overlap = atoi(argv[6]);

	if(overlap != 0 and overlap != 1)
	{
		printf("Overlap is true (1) or false (0)\n");
		exit(1);
	}

	// Vector sequence - Complex numbers
	vector<complex<float> > seq;

	char nuc;	// Input nucleotide, push in sequence
	// Vector of sequence nucleotides
	vector<char> ch_seq;

	// Read from file
	ifstream fin;
	fin.open(argv[1], ios::in);
	// Read first line of the file
	string line;
	getline(fin, line);
	// cout<<"Reading sequence %s\n";
	printf("Reading sequence %s\n", argv[1]);
	// getline(fin, line); // Clear the buffer

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> gen_R(0, 1);
    std::uniform_int_distribution<> gen_Y(0, 1);
    std::uniform_int_distribution<> gen_K(0, 1);
    std::uniform_int_distribution<> gen_M(0, 1);
    std::uniform_int_distribution<> gen_S(0, 1);
    std::uniform_int_distribution<> gen_W(0, 1);
    std::uniform_int_distribution<> gen_B(0, 2);
    std::uniform_int_distribution<> gen_D(0, 2);
    std::uniform_int_distribution<> gen_H(0, 2);
    std::uniform_int_distribution<> gen_V(0, 2);
    std::uniform_int_distribution<> gen_N(0, 3);


	while(fin>>nuc)
	{
		if(nuc != 'A' && nuc != 'T' && nuc != 'C' && nuc != 'G')
		{
			if(nuc == 'R')
			{
				int cur = gen_R(gen);
				if(cur == 0)			
					ch_seq.push_back('A');
				else if(cur == 1)	
					ch_seq.push_back('G');
			}
			else if(nuc == 'Y')
			{
				int cur = gen_Y(gen);
				if(cur == 0)
					ch_seq.push_back('C');
				else if(cur == 1)
					ch_seq.push_back('T');
			}
			else if(nuc == 'K')
			{
				int cur = gen_K(gen);
				if(cur == 0)
					ch_seq.push_back('G');
				else if(cur == 1)
					ch_seq.push_back('T');
			}
			else if(nuc == 'M')
			{
				int cur = gen_M(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('C');
			}
			else if(nuc == 'S')
			{
				int cur = gen_S(gen);
				if(cur == 0)
					ch_seq.push_back('C');
				else if(cur == 1)
					ch_seq.push_back('G');
			}
			else if(nuc == 'W')
			{
				int cur = gen_W(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('T');
				else if(cur == 2)
					ch_seq.push_back('U');
			}
			else if(nuc == 'B')
			{
				int cur = gen_B(gen);
				if(cur == 0)
					ch_seq.push_back('C');
				else if(cur == 1)
					ch_seq.push_back('G');
				else if(cur == 2)
					ch_seq.push_back('T');
			}
			else if(nuc == 'D')
			{
				int cur = gen_D(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('G');
				else if(cur == 2)
					ch_seq.push_back('T');
			}
			else if(nuc == 'H')
			{
				int cur = gen_H(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('C');
				else if(cur == 2)
					ch_seq.push_back('T');
			}
			else if(nuc == 'V')
			{
				int cur = gen_V(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('C');
				else if(cur == 2)
					ch_seq.push_back('G');
			}
			else if(nuc == 'N')
			{
				int cur = gen_N(gen);
				if(cur == 0)
					ch_seq.push_back('A');
				else if(cur == 1)
					ch_seq.push_back('T');
				else if(cur == 2)
					ch_seq.push_back('G');
				else if(cur == 3)
					ch_seq.push_back('C');
			}
		}
		else
			ch_seq.push_back(nuc);
	}
	fin.close();

	printf("Size of input sequence = %ld\n", ch_seq.size());
	// printf("Sequence\n");
	// int rep = 0;
	// for (int i = 0; i < ch_seq.size(); ++i)
	// {		
	// 	if(ch_seq[i] != 'A' && ch_seq[i] != 'T' && ch_seq[i] !='G' && ch_seq[i] != 'C')
	// 		printf("\nOOPS - %c\n", ch_seq[i]);
	// 	else
	// 		printf("%c", ch_seq[i]);
	// 	rep++;
	// 	if(rep == 70)
	// 	{
	// 		rep = 0;
	// 		printf("\n");
	// 	}
	// }
	// exit(0);

	
	// Get the sequence
	if(overlap == 1)
	{
		printf("Overlapping kmers - %d\n", kmer);
		seq = kmer_transform(kmer, ch_seq);
	}
	else
	{
		printf("Disjoint kmers - %d\n", kmer);
		seq = kmer_transform_non_overlap(kmer, ch_seq);
	}
	
	// while(fin>>nuc)
	// {
	// 	if(nuc != 'A' && nuc != 'G' && nuc != 'C' && nuc != 'T')
	// 		continue;
	// 	// cout<<nuc;
	// 	if(nuc == 'A')
	// 	{
	// 		complex<float> cnum(0, 1);
	// 		seq.push_back(cnum);
	// 	}
	// 	else if(nuc == 'G')
	// 	{
	// 		complex<float> cnum(1, 0);
	// 		seq.push_back(cnum);
	// 	}
	// 	else if(nuc == 'T')
	// 	{
	// 		complex<float> cnum(0, -1);
	// 		seq.push_back(cnum);
	// 	}
	// 	else if(nuc == 'C')
	// 	{
	// 		complex<float> cnum(-1, 0);
	// 		seq.push_back(cnum);
	// 	}
	// }	

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
		if(k % 20 == 0)
			printf("Step %d of %d\r", k, win-1);
		fflush(stdout);
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