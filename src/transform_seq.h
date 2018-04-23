#include <iostream>
#include <fstream>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <cmath>
#include <complex>
#include <map>

#define _USE_MATH_DEFINES

// convert from base 10 to base 4, store in k spaces
vector<int> convert_base_k(int num, int k)
{
	int N = num;
	int rem;
	int fac = 1; // Factor rises by 10
	int acc = 0; // Accumulate
	while(N > 0)
	{
		rem = N % 4;
		acc += (fac * rem);
		fac *= 10;
		N /= 4;
	}
	
	vector<int> knum(k, 0);
	int i = 0; int x = acc;
	while(x > 0)
	{
		int rem = x % 10;
		knum[i] = rem;
		x /= 10;
		++i;
	}
	return knum;
}

// Make kmer mapping
map<vector<char>, complex<float> > make_kmap(int k)
{
	vector<char> nucl; // Nucleotides
	nucl.push_back(65);	// A
	nucl.push_back(84);	// T
	nucl.push_back(71);	// G
	nucl.push_back(67);	// C

	// 4^k possibilities - k MUST be SMALL!
	int num_permute = std::pow(4, k);
	vector<vector<char> > permute(num_permute);

	// Get num_permute roots of unity
	float theta = M_PI*2/num_permute;
	// Temporary base 4 number
	vector<int> knum(k);
	// Now create mapping
	map<vector<char>, complex<float> > kmap;

	// Idea - convert to base 4 numbers
	for (int i = 0; i < num_permute; ++i)
	{
		knum = convert_base_k(i, k);
		permute[i].resize(k);
		for (int j = 0; j < k; ++j)
			permute[i][j] = nucl[knum[j]];
		complex<float> cnum(cos(i*theta), sin(i*theta));
		kmap[permutep[i]] = cnum;
	}

	return kmap;
}

// Read in file, and return transformed sequence - kmer change
vector<complex<float> > kmer_transform(int k, vector<char> ch_seq)
{
	int N = ch_seq.size();

	// Seqeunce of complex floats
	vector<complex<float> > seq;
	// Hold kmer
	vector<char> char_map(k);
	// Map to complex value
	map<vector<char>, complex<float> > kmap = make_kmap(k);

	for (int i = 0; i < N - k + 1; ++i)
	{
		for (int j = 0; j < k; ++j)
			char_map[j] = ch_seq[i+j];
		// Get the complex number corresponding to this char map
		seq.push_back(kmap[char_map]);
	}

	return seq;
}