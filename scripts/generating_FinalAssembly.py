#! /usr/bin/env python
#coding: utf-8

##################################################################################################################
#This script is designed for concatenating the scaffolds matched to a particular chromosome (after mummer run)   #
#and generate the chromosome. This script must be applied after running combining_orderedScaffolds.py            #
#Authors: David Peris UW-Madison, Dept Genetics                                                                  #
#                                                                                                                #
##################################################################################################################

import argparse
import pandas as pd #this is how I usually import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import os,sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="folder where outputs of reorderContigsUsingNucmer.R are saved, default = ./", type = str, default = "./")
parser.add_argument("-l","--length", help="size of the scaffold to be accepted, default = 10000", type = int, default = 10000)

args = parser.parse_args()

input_folder = args.input
lenght_scaffold = args.length

strain_folders = glob.glob(input_folder+"*")

#We will rename the reorderContigsUsingNucmer output csv file if the file has been not previously renamed
for strain in strain_folders:
	strain = strain.split("/")[-1]
	if "." not in strain:	#Check the listed file is a folder
		current_folder = input_folder+strain+"/"
		strain_report = open(current_folder+strain+'_reportAssembly.txt', 'w')
		csv_file = pd.read_csv(current_folder+strain+'.csv')
		csv_file = csv_file.dropna(how='all')
		csv_file = csv_file.sort_values(by =['status','ref.start'])
		list_chromosomes = []
		for chromosome_name in csv_file.drop_duplicates(subset=['status'])['status']:
			globals()[chromosome_name] = []
			list_chromosomes.append(chromosome_name)
		list_chromosomes = [x for x in list_chromosomes if str(x) != 'nan']
		list_chromosomes = [x for x in list_chromosomes if str(x) != 'unplaced']
		#We will include the scaffolds to their corresponding chromosome
		for scaffold in csv_file.drop_duplicates(subset=['contig.name'])['contig.name']:
			chromosome_location = csv_file.loc[csv_file['contig.name'] == scaffold,'status'].iloc[0]
			if csv_file.loc[csv_file['contig.name'] == scaffold,'contig.full.length'].iloc[0] > lenght_scaffold:
				globals()[chromosome_location].append(scaffold)
		records = list(SeqIO.parse(current_folder+strain+'_allScaffoldsReordered.fasta','fasta'))
		name_chromosomes = []
		scaffolds_added = []
		for chromosome_name in list_chromosomes:
			strain_report.write("================="+chromosome_name+"=================\n")
			if 'mtDNA' not in chromosome_name.split('_')[-1]:
				name_chr = strain + '_' + chromosome_name.split('_')[-1]
				name_chromosomes.append(name_chr)
				globals()[name_chr] = Seq("", generic_dna)
				for scaffold in globals()[chromosome_name]:
					for sequence in records:
						if scaffold == sequence.id:
							globals()[name_chr] += sequence.seq + ('N' * 5000)
							scaffolds_added.append(scaffold)
							strain_report.write(scaffold+"\n")
		for chromosome_name in name_chromosomes:
			globals()[chromosome_name] = SeqRecord(globals()[chromosome_name][:-5000], id = chromosome_name)
		name_chromosomes = sorted(name_chromosomes)
		sorted_chromosomes = name_chromosomes[:4]+name_chromosomes[5:9]
		sorted_chromosomes.append(name_chromosomes[4])
		sorted_chromosomes += name_chromosomes[9:]
		new_fasta = open(current_folder+strain+'_finalAssembly.fas', 'w')
		new_Allfasta = open(current_folder+strain+'_finalAssembly_wScaffolds.fas', 'w')
		for records in sorted_chromosomes:
			new_fasta.write('>'+globals()[records].id+'\n')
			new_fasta.write(str(globals()[records].seq)+'\n')
			new_Allfasta.write('>'+globals()[records].id+'\n')
			new_Allfasta.write(str(globals()[records].seq)+'\n')
		scaffold_list_noadded = []
		new_fasta.close()
		for scaffold_noadded in SeqIO.parse(current_folder+strain+'_allScaffoldsReordered.fasta','fasta'):
			if scaffold_noadded.id not in scaffolds_added:
				if len(scaffold_noadded.seq) > 5000:
					scaffold_list_noadded.append(scaffold_noadded)
		SeqIO.write(scaffold_list_noadded,new_Allfasta,'fasta')
		new_Allfasta.close()
		strain_report.close()
		print '%s final Assembly done' % (strain)

print 'DONE!'