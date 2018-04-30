# Usage: ./iCorr <input_seq> <window_size> <corr_plot> <slide_dist>

# ./iCorr 

time ./iCorr ../Sequences/NC_000909.fasta 50000 kmer_plots/NC_000909/unamer_50k_10k 10000
time ./iCorr ../Sequences/NC_000909.fasta 10000 kmer_plots/NC_000909/unamer_10k_2k 2000

time ./iCorr ../Sequences/NC_001879.fasta 50000 kmer_plots/NC_001879/unamer_50k_10k 10000
time ./iCorr ../Sequences/NC_001879.fasta 10000 kmer_plots/NC_001879/unamer_10k_2k 2000
