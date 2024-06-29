#!/usr/bin/python3

#======================================================================#
#
#Produce recombination/similarity plots for multiple sequences in a
#fasta alignment
#
#Last Updated: 29 May 2022. anthony.chamings[at]dpi.nsw.gov.au
#
#======================================================================#
import os
import shutil
import sys
import random
import csv
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

window=200
step=20
output_dir=""


def string_similarity(s1,s2):
	score=0
	nt_length=0
	if len(s1)!=len(s2):
		print("Error! Strings must be same size for similarity comparison.")
		print("s1: "+str(s1))
		print("s2: "+str(s2))
		exit(1)
	for pos,value in enumerate(s1):
		'''if s1[pos]=="-" and s2[pos]=="-":
			#this is a blank value
			continue
		if s1[pos]=="-" and s2[pos]!="-":
			nt_length=nt_length+1
		if s1[pos]!="-" and s2[pos]=="-":
			nt_length=nt_length+1'''
		if s1[pos]==s2[pos]:
			if s1[pos]=="N" or s2[pos]=="N":
				score=score+1
			else:
				score=score+1
			nt_length=nt_length+1
		if s1[pos]!=s2[pos]:
		#if s1[pos]!="-" and s2[pos]!="-" and s1[pos]!=s2[pos]:
			nt_length=nt_length+1
	percentage_similarity=round(int(score)/nt_length*100,2)
		
	return percentage_similarity
		
def contains_gaps(s1,s2):
	#Report if string s1 or s2 have any N's present
	s1_count=s1.count("N")
	s2_count=s2.count("N")
	gap_count=s1_count+s2_count
	
	return gap_count

def random_string(rdm_strlen):
	random_string=""
	for i in range(rdm_strlen):
		rdm_chars="abcdefghijklmnopqrstuvwxyz0123456789"
		random_string=random_string+str(rdm_chars[random.randint(0,len(rdm_chars)-1)])
	return random_string

def show_usage():
	print("Usage: recombination_analysis.py -i <fasta alignment> -w <window size> -n <nt step>\n")
	print("\t-i: An alignment in fasta format")
	print("\t-w: size of sliding window in nucleotides. Default: ("+str(window)+")")
	print("\t-n: size of step in nucleoditdes. Default: ("+str(step)+")")
	print("\t-o: output directory to save plots")
	exit(0)

def main(argv):
	global window
	global step
	global output_dir
	arg_count=0
	args=len(argv)
	print(args)
	if "-h" in argv or "-?" in argv or "--help" in argv:
		print("help flag")
		show_usage()
		exit(0)
		
	if args==0:
		show_usage()
		exit(0)
	
	
	
	while arg_count<args:
		arg=argv[arg_count]
		
		if argv[arg_count]=="-i" or argv[arg_count]=="--input":
			alignment_file=argv[arg_count+1]
			arg_count=arg_count+1
		
		if argv[arg_count]=="-w" or argv[arg_count]=="--window":
			window=int(argv[arg_count+1])
			arg_count=arg_count+1
		
		if argv[arg_count]=="-n" or argv[arg_count]=="--step":
			step=int(argv[arg_count+1])
			arg_count=arg_count+1
		if argv[arg_count]=="-o" or argv[arg_count]=="--output":
			output_dir=str(argv[arg_count+1])
			arg_count=arg_count+1
		
		arg_count=arg_count+1
		
	
	#Check if alignment file exists
	if alignment_file=="" or not os.path.isfile(alignment_file):
		print("Cannot find alignment file: "+str(alignment_file))
		exit(1)
	
	#Check if an output_dir is set
	if output_dir =="":
		print("Output directory not specified. Generating random directory.")
		output_dir=str(alignment_file)+"."+str(random_string(8))
		while os.path.isdir(output_dir):
			output_dir=str(alignment_file)+"."+str(random_string(8))

		
	
	
	
	print("Alignment file: "+ str(alignment_file))
	print("Window size: "+str(window))
	print("Step size: "+str(step))
	print("File output will be saved in: "+str(output_dir))




	#delete current output directory and create a new one
	if os.path.isdir(output_dir):
		shutil.rmtree(output_dir)
		#os.rmdir(output_dir)
	os.mkdir(output_dir)





	#read alignment file into memory
	sequences=[]
	sequence_names=[]
	fasta_file=open(alignment_file,"r")
	sequence=""
	sequence_name=""
	for fline in fasta_file:
		if fline[0]==">":
			if sequence_name!="":
				#A new sequence has been reached
				sequence_names.append(sequence_name)
				sequences.append(sequence)
				sequence=""
			sequence_name=fline[1:]
		else:
			sequence=sequence+fline
			
	sequences.append(sequence.upper())
	sequence_names.append(sequence_name)
	fasta_file.close()
	
	#loop through each sequence and compare each sequence
	
	num_sequences=len(sequences)
	
	for s in range(num_sequences):
		current_sequence_name=sequence_names[s].strip()
		current_sequence=sequences[s]
		print("Processing: "+str(current_sequence_name)+"...")
		header=[]
		header.append(current_sequence_name)
		similarity_array=[]
		#gap_array=np.array([])		
		
		#loop through the sequences and compare to the other sequences
		
		for nt in range(0,(len(current_sequence)-window),step):
			
			step_array=[]
			sequence_gap_array=[]
			#because sequence string is zero indexed, add 1 to step
			step_array.append(nt+1)
			sequence_gap_array.append(nt+1)
			
			
			
			
			current_window=current_sequence[nt:(nt+window)]
			for t in range(num_sequences):
				if s!=t:
					target_sequence_name=sequence_names[t].strip()
					target_sequence=sequences[t]
					if (nt==0):
						header.append(target_sequence_name)
						
					target_window=target_sequence[nt:(nt+window)]
					step_array.append(string_similarity(current_window,target_window))
					#TODO - count the gaps for each window
					sequence_gap_array.append(contains_gaps(current_window,target_window))
			similarity_array.append(step_array)
			#gap_array.append(sequence_gap_array)
		#Now put the header at the start of the array
		similarity_array.insert(0,header)
		#gap_array.insert(0,header)
		
		
		sequence_with_highest_similarity=""
		
		outfile=str(output_dir)+"/similarity_"+str(current_sequence_name)+".csv"
		with open(outfile, 'w') as f:
			writer = csv.writer(f, delimiter=',')
			writer.writerows(similarity_array)  
		f.close()
		similarity_df=pd.read_csv(outfile,header=0,index_col=0)
		
		
		# Work out the sequence with the highest similarity and store
		# it in a temporary variable. When that sequence changes, 
		# output the new highest sequence name
		
		
		
		most_similar_sequence_index=similarity_df.idxmax(axis=1)
		most_similar_sequence_value=similarity_df.max(axis=1)
		#print(most_similar_sequence_index,most_similar_sequence_value)
		#similarity_df.to_csv(outfile,sep=",")
		
		
		
		
		
		fig=plt.figure(figsize=(35,8))
		plt.set_cmap("tab20")
		#print("Plotting: "+str(current_sequence_name)+"...")
		
		
		#Prepare the gap array by converting to a percentage which will
		#be used to produce the shade of grey in the background of the 
		#plot
		#print(gap_array)
		#print(gap_array[1:,1:])
		#gap_df=pd.DataFrame(data=gap_array[1:,1:], \
		#					index=gap_array[0,1:], \
		#					columns=gap_array[1:,0])
		
		
		
		for c in similarity_df.columns:
			
			if c in most_similar_sequence_index.values:
				# only plot the most similar sequences for clarity
				# to print all, enable the else statement
				
				line=plt.plot(similarity_df.index.astype(str),similarity_df[c].astype(float),linestyle="-",linewidth=3)
				line[0].set_label(c)
			
			
			else:
				line=plt.plot(similarity_df.index.astype(str),similarity_df[c].astype(float),linestyle="-",linewidth=1)
			
			
			
			
		#plt.legend(labels=similarity_df.columns)
		plt.title(str(current_sequence_name)+" similarity plot",fontsize=32)
		x_lim_min,x_lim_max=plt.xlim()
		
		
		
		
		
		
		plt.xticks(ticks=np.arange(int(x_lim_min), \
			int(x_lim_max), \
			step=1), \
			rotation = 90)
		#plt.xticks(ticks=np.arange(x_lim_min,int(x_lim_max),int(100/step)),rotation = 90)
		
		
		plt.xlabel("Nucleotide position",fontsize=24)
		plt.ylabel("Similarity %",fontsize=24)
		plot_footnote_text="Window size: "+str(window)+"nt , step: "+str(step)+"nt"
		plt.figtext(0.5, -0.05, plot_footnote_text, ha="center", fontsize=18)
		plt.legend(fontsize=18)
		
		
		plt.tight_layout()
		plt.savefig(str(output_dir)+"/similarity_"+str(current_sequence_name)+".png",dpi=150,bbox_inches='tight')
		plt.close()
		print("Done.")
	
	
	



if __name__ == "__main__":
	argv=sys.argv[1:]								#get the arguments passed to this script
	main(argv)
	
