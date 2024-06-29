# similarity_plots
A python script for creating similarity plots of sequences within a FASTA alignment. 

This scripts relies on a number of libraries, some of which are not installedby default when Python is first installed:

Default libraries:
os, shutil, sys, random, csv

Installed libaries:
pandas, matplotlib, numpy 

This script will output a plot for each sequence comparing it against all others within the aligned FASTA file. Recombinant sequences can be detected by looking for sequences which have sudden changes in percentage similarity between one or more aligned sequences. 

The default setting is a 200nt window with 20nt steps but these can be altered viaarguments when the script is run.

Usage: recombination_analysis.py -i <fasta alignment> -w <window size> -n <nt step>\n -o <output directory>
	-i: An alignment in fasta format
	-w: size of sliding window in nucleotides. Default: 200
	-n: size of step in nucleoditdes. Default: 20
	-o: output directory in which plots are saved
