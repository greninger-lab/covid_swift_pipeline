# This file takes the bcftools_variants.csv and fixes the complex mutations.

import subprocess 
import argparse
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd
import sys
import math

def correct_deletions(df, group_to_correct, correct_depth):
	for row_num, row in (group_to_correct[group_to_correct['TCOV'] <= 0.4* correct_depth]).iterrows():
		df.at[row_num,'TCOV'] = correct_depth
		#df.at[row_num,'AAFREQ'] = row['VCOV'] / correct_depth
	return df 

def translate(seq):
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 
		'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W', 
	} 
	protein = "" 
	if "-" in seq:
		return "fs"
	if len(seq)%3 == 0: 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			protein+= table[codon] 
	return protein 

def process_fasta(fasta):
	gene_info = pd.DataFrame(columns=['gene_name','gene_start','nt_seq','aa_seq'])
	gene_name = ""
	gene_start = 0
	for line_num,line in enumerate(open(fasta)):	
		if line_num % 2 == 0:
			gene_name = line.split("transcript:")[1].split("Comment")[0].rstrip()
			gene_start = int(line.split(")")[0].split(":")[-1])
		else:
			gene_info.loc[line_num-1] = [gene_name] + [gene_start] + [line.rstrip()] + [translate(line.rstrip())]
	return gene_info

def find_new_nts(aa_nts, list_of_snpids, first_nt_pos):
	new_nts = list(aa_nts)
	for snpid in list_of_snpids:
		nt_pos = int(snpid[1:-1])
		new_nts[nt_pos - first_nt_pos] = snpid[-1]
	return(''.join(new_nts))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-name', help="Provide sample name.")
	parser.add_argument('-file', help="Provide processed bcftools_variants.csv file.")
	parser.add_argument('-fasta', help="Provide AT_refGeneMrna.fa file with gene sequences.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	sample_name = args.name
	input_variants = args.file

	gene_info = process_fasta(args.fasta)

	output_variants_file = sample_name + "_bcftools_variants_edited.csv"

	df = pd.read_csv(input_variants).fillna("-")

	# First, let's correct some variants that are obviously wrong
	# Grabbing variants at the same aa position
	for name,group in df.groupby(["gene","AAPOS"]):
		if group.shape[0] >1:
			for nuc,nuc_group in group.groupby(["NTPOS"]):
				# Here we are protecting against real deletions, which mess up depths/afs of surrounding variants
				correct_depth = max(nuc_group["TCOV"])
				df = correct_deletions(df, nuc_group, correct_depth)
				# print(df.iloc[row_num])
			
			# Here we protect against whole amino acid deletions, such as the spike 143-144 deletion
			if "-" in group['AASUB'].unique():
				max_depth = max(group['TCOV'])
				real_depth = group[group['AASUB']=="-"]['TCOV'].values[0]
				if real_depth == max_depth:
					df = correct_deletions(df, group, real_depth)

	# Now we recalculate AAFREQ and refilter by 0.01
	df['AAFREQ'] = df['VCOV'] / df['TCOV']
	df = df.loc[df['AAFREQ'] >= 0.01]

	# Now for the actual complex mutations part!
	for name,group in df.groupby(["gene","AAPOS"]):
		# We don't want to bother with super low-frequency variants
		filtered_group = group.loc[group['AAFREQ'] >= 0.05]
		# Throw out those with deletions
		filtered_group = group.loc[group['AASUB']!= "-"]
		if filtered_group.shape[0] > 1:
			filtered_group = filtered_group.sort_values(['AAFREQ'], ascending=False)
			for i in range(len(filtered_group.index) - 1):
				j = i+1
				list_to_combine = pd.DataFrame()
				list_to_combine = list_to_combine.append(filtered_group.iloc[i])
				given_nuc = filtered_group.iloc[i]['NTPOS']

				# This list is sorted from high to low AAFREQ. So go through and add mutations that could combine with
				# the given mutation at i that are within 5% frequency.
				while(j < len(filtered_group.index) and filtered_group.iloc[i]['AAFREQ'] - filtered_group.iloc[j]['AAFREQ'] <=0.05):
					compare_nuc = filtered_group.iloc[j]['NTPOS']
					if given_nuc != compare_nuc:
						list_to_combine = list_to_combine.append(filtered_group.iloc[j])

					j+=1

				# If there are multiple mutations at the same nucleotide, we have no way
				# of telling which one is the complex one. So we drop these and report separately
				list_to_combine = list_to_combine.drop_duplicates(subset=['NTPOS'],keep=False).sort_values(['NTPOS'])

				# Ok we're done filtering... now the actual complex mutation part
				if list_to_combine.shape[0] > 1:
					print(list_to_combine)
					# Grab relevant gene info
					rel_gene = gene_info.loc[gene_info['gene_name'] == name[0]]

					# Grab the nt position of first nt in relevant aa
					nt_pos_in_gene = int(name[1]) * 3
					first_nt_pos = int(nt_pos_in_gene) + int((rel_gene['gene_start'] - 2).values[0])
					if name[0] == "ORF1ab_polyprotein_ribosomal_slippage":
						first_nt_pos = first_nt_pos - 13203
						nt_pos_in_gene = nt_pos_in_gene - 13203

					# Grab nts
					aa_nts = rel_gene['nt_seq'].values[0][nt_pos_in_gene-3:nt_pos_in_gene]
					
					# Grab new changes
					new_nts = find_new_nts(aa_nts,(list(list_to_combine['snpid'])), first_nt_pos)
					new_aa = translate(new_nts)
					new_snpid = aa_nts + str(first_nt_pos) + new_nts 
					
					# Implement new changes in original dataframe
					row_to_change = list_to_combine.iloc[0].name
					df.at[row_to_change,'NSPSUB']=new_aa
					df.at[row_to_change,'AASUB']=new_aa
					df.at[row_to_change,'snpid']=new_snpid

					# Get rid of other rows old dataframe
					for i in range(1,list_to_combine.shape[0]):
						df = df.drop([list_to_combine.iloc[i].name], errors="ignore")
					
					#print(df.loc[df['AAPOS']==name[1]])

	#print(df[df['AAPOS']==144])

	# Drop rows with empty values - happens with 3-nucleotide complex changes that have already been covered by the 2-nt case
	df = df.dropna(subset=['Sample'])

	df.to_csv(output_variants_file, index=False)
