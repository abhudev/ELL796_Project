#include <iostream>
#include <fstream>
// #include <bits/stdc++.h>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

typedef long long ll;

int main(int argc, char const *argv[])
{
	if(argc != 5)
	{
		printf("Usage: ./GC <input_seq> <window_size> <GC_plot> <slide>\n");
		exit(1);
	}

	// Window
	int win = atoi(argv[2]);
	int slide = atoi(argv[4]);

	// Vector sequence - G is 1, C is 
	vector<char> seq;

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
		seq.push_back(nuc);
	}

	fin.close();

	int N = seq.size();
	cout<<"Length of sequence is "<<N<<endl;

	// Compute number of windows - considering only 'whole' windows
	int num_win = N - win + 1;

	// Array to hold the cumulative GC skew values
	vector<float> cum_gc(num_win);

	// Now compute. First get for the first window, then keep updating
	float G = 0, C = 0;
	int i;
	for(i = 0; i < win; ++i)
	{
		switch(seq[i])
		{
			case 'G':
						++G;
						break;
			case 'C':	
						++C;
						break;
			default:
						break;
		}
	}

	// Update first value of cumulative GC skew
	cum_gc[0] = (C-G)/(G+C);
	cout<<endl;
	// Now compute for all other windows
	for(i = 1; i < num_win; ++i)
	{
		// if((i % 100) == 0)
			// cout<<"\rStep "<<i;
		// switch(seq[i - 1])
		// {
		// 	case 'G':
		// 			 --G;
		// 			 break;
		// 	case 'C':
		// 			--C;
		// 			break;
		// 	default:
		// 			break;
		// }
		if(seq[i-1] == 'G')
			--G;
		else if(seq[i-1] == 'C')
			--C;
		// switch(seq[i + win - 1])
		// {
		// 	case 'G':
		// 				++G;
		// 				break;
		// 	case 'C':	
		// 				++C;
		// 				break;
		// 	default:
		// 				break;
		// }
		if(seq[i+win-1] == 'G')
			++G;
		else if(seq[i+win-1] == 'C')
			++C;
		// Update cumulative skew
		if(G-C != 0)
			cum_gc[i] =  (C-G)/(G+C);
		else
			cum_gc[i] = 0;
	}
	cout<<"\nDone\n";
	// Write cumulative GC to file
	FILE* fout;
	fout = fopen(argv[3], "w");
	// Increment according to slide 
	for(i = 0; i < num_win; i += slide)
	{		
		fprintf(fout, "%f\n", cum_gc[i]);
	}
	cout<<endl;
	fclose(fout);
	return 0;
}