

"""
Author: Ruth Hanna
Email: rhanna@broadinstitute.org
Title: Identifying PAMs near ClinVar SNPs
Description: This script reads in two input files:
1. A raw file of SNP positions downloaded from the ClinVar website.
2. A file with a list of Cas nucleases and information about their PAMs and editing windows.
The script first identifies the genomic context in the reference genome for each SNP and then searches for PAMs for each Cas in the genomic context.
The output is a list of sgRNAs, annotated with the SNP and Cas, and a list of SNPs, annotated with whether they are editable by each Cas.
All genomic locations use 1-based indexing; all locations in sequences (e.g. PAM_Location) use 0-based indexing.
"""

import pandas as pd
import numpy as np
import argparse, csv, itertools, gzip, os, re
from datetime import datetime
import matplotlib.pyplot as plt
plt.switch_backend('agg')


BUFFER_LEN = 50


nucleotide_codes_dict = {
	'A':{'A'},
	'C':{'C'},
	'G':{'G'},
	'T':{'T'},
	'R':{'A','G'},
	'Y':{'C','T'},
	'W':{'A','T'},
	'S':{'G','C'},
	'M':{'A','C'},
	'K':{'G','T'},
	'B':{'G','C','T'},
	'H':{'A','C','T'},
	'D':{'A','G','T'},
	'V':{'A','G','C'},
	'N':{'A','G','C','T'}}


complements_dict = {
	'A':'T',
	'T':'A',
	'C':'G',
	'G':'C',
	'N':'N',
	'R':'Y',
	'Y':'R',
	'W':'W',
	'S':'S',
	'M':'K',
	'K':'M',
	'B':'V',
	'V':'B',
	'H':'D',
	'D':'H'}



def GetParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant_file',
        type=str,
        help='.txt input file from the ClinVar website)')
    parser.add_argument('--pam_file',
    	type=str,
    	help='.txt input file with the following columns: name of Cas, PAM, start of editing window, end of editing window, sgRNA length, allowed nt (e.g. for base editing)')
    parser.add_argument('--output_name',
    	type =str,
    	help='name of output file')
    parser.add_argument('--sgrna_type',
    	type=str,
    	choices=['edit_in','edit_out'],
    	help='edit_in for sgRNAs to introduce ClinVar SNPs (i.e. sgRNA targets WT sequence), edit_out for sgRNAs to correct ClinVar SNPs (i.e. sgRNA targets mutant sequence)')
    parser.add_argument('--pathogenic_only',
    	type=str,
    	choices=['True','False'],
    	help='whether or not to only include pathogenic SNPs in ClinVar')
    parser.add_argument('--hdr_distance',
    	type=int,
    	help='distance to allow between cut site and SNP for HDR-based editing')    
    return parser



def ParseVariantFile(variant_df,pathogenic_only):
	# Cleans variant file from ClinVar database.
	# Removes non-GRCh38 rows, non-SNPs, and chromosomes other than 1-22, X, and Y.
	# Returns cleaned df.
	
	print('Parsing variant file')
	if pathogenic_only == 'True':
		parsed_variant_df = variant_df[
			(variant_df.Assembly == 'GRCh38') &
			(variant_df.ClinicalSignificance == 'Pathogenic') & 
			(variant_df.Type == 'single nucleotide variant') & 
			(variant_df.Chromosome != 'MT') & 
			(variant_df.Chromosome != 'na') & 
			(variant_df.ReferenceAllele != 'na')]
	else:
		parsed_variant_df = variant_df[
			(variant_df.Assembly == 'GRCh38') &
			(variant_df.Type == 'single nucleotide variant') & 
			(variant_df.Chromosome != 'MT') & 
			(variant_df.Chromosome != 'na') & 
			(variant_df.ReferenceAllele != 'na')]
	parsed_variant_df.index = range(0,len(parsed_variant_df))
	return parsed_variant_df



def ParsePAMFile(pam_df,hdr_distance):
	# Converts input PAMs from a string to a list.
	split_row_list = []
	column_list = list(pam_df.columns.values)
	column_list[column_list.index('PAM')] = 'Split_PAM'
	for row in pam_df['PAM']:
		split_row = row.split(',')
		split_row_list.append(split_row)
	pam_df['Split_PAM'] = split_row_list
	pam_df = pam_df.drop(['PAM'],axis=1)
	pam_df = pam_df[column_list]
	pam_df['Edit_Window_Start'] = pam_df['Edit_Window_Start'].apply(lambda x: x - hdr_distance)
	pam_df['Edit_Window_End'] = pam_df['Edit_Window_End'].apply(lambda x: x + hdr_distance)
	return pam_df



def ReadChromFiles(start_chrom,end_chrom):
	# Loads raw chromosome files and stores as a dictionary.

	chr_dict = {}
	for chrom_num in range(start_chrom,end_chrom+1):
		if chrom_num == 23:
			chrom_num = 'X'
		if chrom_num == 24:
			chrom_num = 'Y'
		fasta_path = '/rnai/refdb/genome/human/GRCh38/'
		chr_file = fasta_path + 'hs_ref_GRCh38_chr' + str(chrom_num) + '.fa.gz'
		print('Reading chrom ' + str(chrom_num))
		f = gzip.open(chr_file,'r')
		chr_seq = f.read()
		f.close()
		chr_seq_split = chr_seq.split('\n')
		chr_seq = ''.join(chr_seq_split[1:])
		chr_dict[str(chrom_num)] = chr_seq
	return chr_dict



def LocateSNP(chr_dict, parsed_variant_df, sgrna_type):
	# Takes in a df of SNP genomic locations and the dictionary of chromosome sequences.
	# For each SNP, extracts the 50 nts preceding the SNP and the 50 nts following the SNP.
	# Outputs a snp_df with the genomic context of the given SNP.

	snp_location_list = []
	error_list = []
	snp_df_columns_list = list(parsed_variant_df.columns.values)
	snp_df_columns_list.extend(['Genomic_Context_Start','Genomic_Context_End','SNP_Position','Preceding_Genomic_Context','Succeeding_Genomic_Context','Genomic_Context'])
	for index, row in parsed_variant_df.iterrows():
		r = list(row)
		snp_chr_seq = chr_dict[row['Chromosome']]
		snp_position = row['Start']
		# Extract surrounding genomic context (taking into account 0-based indexing in Python).
		preceding_genomic_context = snp_chr_seq[(snp_position-BUFFER_LEN-1):snp_position-1]
		succeeding_genomic_context = snp_chr_seq[(snp_position):(snp_position+BUFFER_LEN)]
		# For edit_in, design sgRNAs against reference allele.
		if sgrna_type == 'edit_in':
			snp_genomic_context = preceding_genomic_context + row['ReferenceAllele'] + succeeding_genomic_context
		# For edit_out, design sgRNAs against alternate allele.
		elif sgrna_type == 'edit_out':
			snp_genomic_context = preceding_genomic_context + row['AlternateAllele'] + succeeding_genomic_context
		
		r.extend([snp_position-BUFFER_LEN, snp_position+BUFFER_LEN, snp_position, preceding_genomic_context, succeeding_genomic_context, snp_genomic_context])
		snp_location_list.append(r)
		# Check that reference alleles match. If they don't, append to error_list.
		if snp_chr_seq[snp_position-1] != row['ReferenceAllele']:
			error_list.append(r)
		# Check that length of extracted genomic sequence is correct. If it isn't, append to error_list.
		if len(snp_genomic_context) != 2*BUFFER_LEN + 1:
			error_list.append(r)
	snp_df = pd.DataFrame(data=snp_location_list,columns = snp_df_columns_list)
	snp_error_df = pd.DataFrame(data=error_list,columns = snp_df_columns_list)
	return snp_error_df, snp_df


	
def ReverseComplement(sequence):
	# Returns the reverse complement of input sequence.

	nucleotides = list(sequence)
	complement = [complements_dict[nucleotide] for nucleotide in nucleotides]
	complement = ''.join(complement)
	reverse_complement = complement[::-1]
	return reverse_complement



def FindGuides(parsed_pam_df, snp_df,sgrna_type):
	# Takes in a parsed_pam_df and a snp_df.
	# Find guides for each Cas in the genomic context of each SNP.
	# Outputs a two-element tuple with a sgrna_df (indexed by sgRNA) and a summary_df (indexed by SNP).

	pam_location_list = []
	summary_list = []
	# Create a list of column headers for parsed_pam_df, FindGuides input values, and sgrna_df.
	snp_df_columns_list = list(snp_df.columns.values)
	pam_df_columns_list = list(parsed_pam_df.columns.values)
	input_values_labels = snp_df_columns_list[:]
	input_values_labels.extend(pam_df_columns_list)
	input_values_labels.append('Strand')
	sgrna_df_columns_list = input_values_labels[:]
	sgrna_df_columns_list.extend(['PAM_Window','sgRNA','PAM_Location','Found_PAM','Editing_Window','Number_Bystander_Edits','Edit Nucleotides'])

	for i in range(0,len(snp_df)):
		r = list(snp_df.iloc[i,:])
		# pam_count_list stores the pam counts for each row.
		pam_count_list = r[:]
		# The pam counts for each Cas are appended to temp_pam_count_list.
		temp_pam_count_list = []
		# Create a temporary dictionary to map position in sequence to position in ref genome.
		sequence = r[snp_df_columns_list.index('Genomic_Context')]
		temp_genomic_location_dict = {}
		for base in range(0,len(sequence)):
			temp_genomic_location_dict[base] = r[snp_df_columns_list.index('Genomic_Context_Start')] + base
		# Go through each Cas in parsed_pam_df.
		for index in range(0,len(parsed_pam_df)):
			row = list(parsed_pam_df.iloc[index,:])
			# Initialize PAM count to 0.
			pam_count = 0
			# Create a list of input values for FindGuides in sense strand.
			input_values = r[:]
			input_values.extend(row)
			# Search for guides in sense strand.
			input_values_sense = input_values[:]
			input_values_sense.append('sense')
			find_guides_output_sense = FindGuidesInSequence(input_values_sense,input_values_labels,temp_genomic_location_dict,sgrna_type)
			for guide in range(0,len(find_guides_output_sense)):
				# If a guide is found, increment the PAM count.
				pam_count += 1
				pam_location_list.append(find_guides_output_sense[guide])
			# Create a list of input values for FindGuides in antisense strand.
			input_values_antisense = input_values[:]
			input_values_antisense.append('antisense')
			input_values_antisense[input_values_labels.index('Genomic_Context')] = ReverseComplement(r[snp_df_columns_list.index('Genomic_Context')])
			input_values_antisense[input_values_labels.index('ReferenceAllele')] = ReverseComplement(r[snp_df_columns_list.index('ReferenceAllele')])
			input_values_antisense[input_values_labels.index('AlternateAllele')] = ReverseComplement(r[snp_df_columns_list.index('AlternateAllele')])
			# Search for guides in antisense strand.
			find_guides_output_antisense = FindGuidesInSequence(input_values_antisense,input_values_labels,temp_genomic_location_dict,sgrna_type)
			for guide in range(0,len(find_guides_output_antisense)):
				# Reset ref allele, alt allele, and genomic context.
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('ReferenceAllele')] = r[snp_df_columns_list.index('ReferenceAllele')]
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('AlternateAllele')] = r[snp_df_columns_list.index('AlternateAllele')]
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('Genomic_Context')] = r[snp_df_columns_list.index('Genomic_Context')]
				# If a guide is found, reset PAM_Location (still indicating the first nucleotide in the PAM) and increment pam_count.
				find_guides_output_antisense[guide][sgrna_df_columns_list.index('PAM_Location')] = len(sequence)-1-find_guides_output_antisense[guide][sgrna_df_columns_list.index('PAM_Location')]
				pam_count += 1
				pam_location_list.append(find_guides_output_antisense[guide])
			temp_pam_count_list.append(pam_count)
		# Get total number of PAM sites.
		pam_count_list.extend(temp_pam_count_list + [sum(temp_pam_count_list)])
		summary_list.append(pam_count_list)
	# Create summary dataframe indexed by SNP.
	summary_df_columns_list = snp_df_columns_list[:]
	summary_df_columns_list.extend(parsed_pam_df['Cas'])
	summary_df_columns_list.extend(['PAM_Count'])	
	if summary_list:
		summary_df = pd.DataFrame(data=summary_list, columns=summary_df_columns_list)
	if not summary_list:
		print('Error: summary_list is empty')
		return
	# Create sgRNA dataframe indexed by sgRNA.
	if pam_location_list:
		sgrna_df = pd.DataFrame(data=pam_location_list, columns=sgrna_df_columns_list)
	if not pam_location_list:
		print('Error: pam_location_list is empty')
		return
	return sgrna_df, summary_df



def FindGuidesInSequence(input_values, input_values_labels, temp_genomic_location_dict, sgrna_type):
	# Takes in information about SNP, sequence, and PAM.
	# Finds guides using CountPAMs.
	# Outputs a list of guides in given sequence.

	sgrna_list = []
	allowed_ref_allele,allowed_alt_allele,allowed_ref_allele_set,allowed_alt_allele_set,chrom_num,spacer_len,sequence,edit_window_start,edit_window_end,pam_list,pamside,strand,ref_allele,alt_allele,snp_position,genomic_context_start = GetInputValues(input_values,input_values_labels,sgrna_type)
	# Check whether nucleotide is allowed.
	if ref_allele not in allowed_ref_allele_set:
		return sgrna_list
	if alt_allele not in allowed_alt_allele_set:
		return sgrna_list
	# Calculate snp_position_in_sequence (should equal BUFFER_LEN).
	snp_position_in_sequence = snp_position - genomic_context_start
	if snp_position_in_sequence != BUFFER_LEN:
		print('Error: snp_position_in_sequence != BUFFER_LEN')
	# Go through each PAM in the list of input PAMs.
	for pam in pam_list:
		pam_window, pam_window_start, pam_window_end = GetPamWindow(edit_window_start,edit_window_end,spacer_len,pamside,snp_position_in_sequence,sequence,pam)
		# Check whether given genomic context is long enough for pam window.
		if snp_position_in_sequence+pam_window_start < 0 or snp_position_in_sequence+pam_window_end+len(pam) > len(sequence)-1:
			print('Error: increase BUFFER_LEN')
			return sgrna_list
		# Search for pam in strand.
		pams_in_pam_window_list = CountPAMs(pam_window,pam)
		# Extract guides.
		sgrna_list = ExtractSpacers(pams_in_pam_window_list,sgrna_list,input_values,snp_position_in_sequence,pam_window,pam_window_start,pamside,sequence,spacer_len,pam,edit_window_start,edit_window_end,sgrna_type,allowed_alt_allele,allowed_ref_allele)
				
	return sgrna_list



def ExtractSpacers(pams_in_pam_window_list,sgrna_list,input_values,snp_position_in_sequence,pam_window,pam_window_start,pamside,sequence,spacer_len,pam,edit_window_start,edit_window_end,sgrna_type,allowed_alt_allele,allowed_ref_allele):
	if pams_in_pam_window_list:
		for pam_location_list in pams_in_pam_window_list:
			num_bystander = 0
			edit_nucs = ''
			pam_location = pam_location_list[0]
			found_pam = pam_location_list[1]
			pam_location_in_sequence = pam_location + snp_position_in_sequence + pam_window_start
			if pamside == 3:
				sgrna = sequence[(pam_location_in_sequence-spacer_len):pam_location_in_sequence]
			elif pamside == 5:
				sgrna = sequence[(pam_location_in_sequence+len(pam)):(pam_location_in_sequence+len(pam)+spacer_len)]
			# Check for non-ACGT characters
			if re.search(r"[^ACGT]"	,sgrna):
				error = 'Found non-ACGT base'
				print(error)
			editing_window = sgrna[(edit_window_start - 1):(edit_window_end)]
			if sgrna_type == 'edit_out':
				if allowed_alt_allele != 'N':
					num_bystander = editing_window.count(allowed_alt_allele) - 1
					edit_nucs = GetEditNucs(sgrna,allowed_alt_allele,edit_window_start)
			elif sgrna_type == 'edit_in':
				if allowed_ref_allele != 'N':
					num_bystander = editing_window.count(allowed_ref_allele) - 1
					edit_nucs = GetEditNucs(editing_window,allowed_ref_allele,edit_window_start)
			sgrna_list.append(input_values + [pam_window,sgrna,pam_location_in_sequence,found_pam,editing_window,num_bystander,edit_nucs])
	return sgrna_list



def GetInputValues(input_values,input_values_labels,sgrna_type):
	edit_from = input_values[input_values_labels.index('Edit_From')]
	edit_to = input_values[input_values_labels.index('Edit_To')]
	if sgrna_type == 'edit_in':
		allowed_ref_allele = edit_from
		allowed_alt_allele = edit_to
	elif sgrna_type == 'edit_out':
		allowed_ref_allele = edit_to
		allowed_alt_allele = edit_from
	allowed_ref_allele_set = nucleotide_codes_dict[allowed_ref_allele]
	allowed_alt_allele_set = nucleotide_codes_dict[allowed_alt_allele]
	chrom_num = input_values[input_values_labels.index('Chromosome')]
	spacer_len = input_values[input_values_labels.index('Spacer_Len')]
	sequence = input_values[input_values_labels.index('Genomic_Context')]
	edit_window_start = input_values[input_values_labels.index('Edit_Window_Start')]
	edit_window_end = input_values[input_values_labels.index('Edit_Window_End')]
	pam_list = input_values[input_values_labels.index('Split_PAM')]
	pamside = input_values[input_values_labels.index('PAMside')]
	strand = input_values[input_values_labels.index('Strand')]
	ref_allele = input_values[input_values_labels.index('ReferenceAllele')]
	alt_allele = input_values[input_values_labels.index('AlternateAllele')]
	snp_position = input_values[input_values_labels.index('SNP_Position')]
	genomic_context_start = input_values[input_values_labels.index('Genomic_Context_Start')]
	return allowed_ref_allele,allowed_alt_allele,allowed_ref_allele_set,allowed_alt_allele_set,chrom_num,spacer_len,sequence,edit_window_start,edit_window_end,pam_list,pamside,strand,ref_allele,alt_allele,snp_position,genomic_context_start



def GetPamWindow(edit_window_start,edit_window_end,spacer_len,pamside,snp_position_in_sequence,sequence,pam):
	# Slice out pam window, taking into account the length of the PAM.
	if pamside == 3:
		pam_window_start = spacer_len - edit_window_end + 1
		pam_window_end = spacer_len - edit_window_start + len(pam) + 1
	if pamside == 5:
		pam_window_start = -1*edit_window_end - len(pam) + 1
		pam_window_end = -1*edit_window_start + 1
	pam_window = sequence[(snp_position_in_sequence+pam_window_start):(snp_position_in_sequence+pam_window_end)]
	return pam_window, pam_window_start, pam_window_end



def GetEditNucs(editing_window, nuc, edit_window_start):
	edit_nucs = ''
	for i,pos in enumerate(editing_window):
		if pos == nuc:
			edit_nucs += nuc+'_'+str(i+edit_window_start)+';'
	return edit_nucs



def CountPAMs(sequence,pam):
	# Takes in sequence (5' to 3') and a PAM to search for.
	# Outputs a list of PAM locations, where the location indicates the 5' end of the PAM.
	# CountPAMs is agnostic to strand.
	# To find PAMs in the antisense strand, CountPAMs takes in the reverse complement (5' to 3') and outputs the 5' location of the PAM in the antisense strand.

	pam_location_list = []
	for base in range(0, len(sequence)+1-len(pam)):
		pam_found = True
		for i,nuc in enumerate(pam):
			# Check whether nucleotide is in the set of allowed nucleotides for the PAM
			if sequence[base+i] not in nucleotide_codes_dict[nuc]:
				pam_found = False
				break
		if pam_found == True:
			pam_location_list.append([base,sequence[base:base+len(pam)]])
	return pam_location_list



def GetStats(summary_df,cas_list,output_folder,output_name):
	os.makedirs(output_folder+'/figures')
	for cas in cas_list:
		fig,ax=plt.subplots()
		plt.hist(summary_df[cas],bins=range(min(summary_df[cas]), max(summary_df[cas])+3, 1))
		plt.xlabel('number of sgRNAs')
		plt.ylabel('number of SNPs')
		plt.title(cas)
		fig.savefig(output_folder + '/figures/' + cas + '.pdf')
		plt.close()
	return



def main():
	# Read in and parse variant file from ClinVar.
	args = GetParser().parse_args()
	variant_df = pd.read_table(args.variant_file, dtype = {'#AlleleID':str,'Name':str,'ClinicalSignificance':str,'Assembly':str,'Chromosome':str,'Type':str,'Start':int,'ReferenceAllele':str,'AlternateAllele':str,'ReviewStatus':str}, usecols = ['#AlleleID','Name','ClinicalSignificance','Assembly','Chromosome','Type','Start','ReferenceAllele','AlternateAllele','ReviewStatus'])
	# Parse variant_df.
	parsed_variant_df = ParseVariantFile(variant_df,args.pathogenic_only)
	print(len(parsed_variant_df))
	# Read in chrom files.
	chr_dict = ReadChromFiles(1,24)
	# Locate SNPs.
	snp_error_df, snp_df = LocateSNP(chr_dict, parsed_variant_df,args.sgrna_type)
	# Read in pam file and parse list of PAMs into list. 
	# Add HDR distance to edit windows for each PAM
	parsed_pam_df = ParsePAMFile(pd.read_table(args.pam_file),args.hdr_distance)

	# Count pams in each snp and produce list of sgRNAs.
	sgrna_df, summary_df = FindGuides(parsed_pam_df,snp_df,args.sgrna_type)
	output_folder = args.output_name + '_' + str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	GetStats(summary_df,parsed_pam_df['Cas'],output_folder,args.output_name)
	# Write files.
	parsed_variant_df.to_csv(output_folder + '/parsed_variant_df_' + args.output_name + '.txt',sep='\t',index=False)
	snp_df.to_csv(output_folder + '/snp_df_' + args.output_name + '.txt',sep='\t',index=False)
	if not snp_error_df.empty:
		snp_error_df.to_csv(output_folder + '/snp_error_df_' + args.output_name + '.txt',sep='\t',index=False)
	sgrna_df.to_csv(output_folder + '/sgrna_df_' + args.output_name + '.txt', sep='\t',index=False)
	summary_df.to_csv(output_folder + '/summary_df_' + args.output_name + '.txt', sep='\t',index=False)
	with open(output_folder+'/README.txt','w') as o:
		w = csv.writer(o,delimiter='\t')
		w.writerow((['Variant file: '+args.variant_file]))
		w.writerow((['PAM file: '+args.pam_file]))
		w.writerow((['Output folder: '+output_folder]))
		w.writerow((['HDR distance: '+str(args.hdr_distance)]))



if __name__ == '__main__':
	main()