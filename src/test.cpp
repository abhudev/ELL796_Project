#include <iostream>
#include <fstream>
#include <vector>
// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <complex>
#include "transform_seq.h"


int main(int argc, char const *argv[])
{
	if(argc != 3)
	{
		printf("Usage: ./test <infile> <kmer>\n");
		exit(1);
	}

	int kmer = atoi(argv[2]);

	ifstream fin;
	fin.open(argv[1], ios::in);
	// Read first line of the file
	string line;
	getline(fin, line);
	cout<<"Reading sequence\n";

	char nuc;	// Input nucleotide, push in sequence
	// Vector of sequence nucleotides
	vector<char> ch_seq;

	vector<complex<float> > seq;

	while(fin>>nuc)
		ch_seq.push_back(nuc);
	fin.close();

	// Get the sequence
	seq = kmer_transform(kmer, ch_seq);

	printf("Sizes of sequences are %ld and %ld\n", ch_seq.size(), seq.size());

	printf("First few items of sequence:\n");
	for (int i = 0; i < 10; ++i)
		printf("%f, %f\n", seq[i].real(), seq[i].imag());

	return 0;
}