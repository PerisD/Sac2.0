#! /usr/bin/env python
#coding: utf-8

##################################################################################################################
#This script is designed for locating reordered and not reordered fasta sequences in one file                    #
#This script must be applied after running reorderContigsUsingNucmer.R                                           #
#Authors: David Peris UW-Madison, Dept Genetics                                                                  #
#Usage: python combining_orderedScaffolds.py ARG1 ARG2                                                           #
#                                                                                                                #
#ARG1: folder where output of reorderContigsUsingNucmer.R are saved                                              #
#ARG2: Path to the original textfile  used in runMummer_serially.py: tabulated text file 1st column              #
#(Path2Assembly/Assembly name) (Path2Reference/Reference_Name)                                                   #
#                                                                                                                #
##################################################################################################################

import sys,os
from Bio import SeqIO
import glob

input_folder = sys.argv[1]
input_file = sys.argv[2]

strain_folders = glob.glob(input_folder+"*") #List the folders contained in the

def CombineReorderedScaffolds(strain): 
	cmdA = 'cat '+current_folder+'*.fas* >' + current_folder+strain + '_reordered.fasta'
	return cmdA

def EMBOSS_info(assembly_path,strain):
	cmdB = 'infoseq ' + assembly_path + ' --auto -only -name -length -pgc > ' + strain + '_infoseq.txt'
	return cmdB

def all_contigs(strain): #Command still not tested
	cmdC = 'cat ' + current_folder + strain + '*eordered.fasta > ' + current_folder+strain+'_allScaffoldsReordered.fasta'
	return cmdC

#We will rename the reorderContigsUsingNucmer output csv file if the file has been not previously renamed
for strain in strain_folders:
	strain = strain.split("/")[-1]
	if "." not in strain:	#Check the listed file is a folder
		current_folder = input_folder+strain+"/"
		report_file = open(strain+'_report_generating_correctAssembly.txt','w')
		#os.rename(current_folder+"remap.coords.result.csv",current_folder+strain+".csv")
		os.system(CombineReorderedScaffolds(strain))
		scaffolds_reorderedFile = []
		for seqRecord_reorderedFile in SeqIO.parse(current_folder+strain + '_reordered.fasta',"fasta"):
			scaffolds_reorderedFile.append(seqRecord_reorderedFile.name)
		scaffolds_originalAssembly = []
		list_of_assemblies = open(input_file)
		for line in list_of_assemblies:
			line = line.strip()
			assembly_line = line.split('\t')[0].split('/')[-1]
			if strain in assembly_line:
				for seqRecord_originalAssembly in SeqIO.parse(line.split('\t')[0],"fasta"):
					scaffolds_originalAssembly.append(seqRecord_originalAssembly)
				report_file.write('============================ %s ====================================\n' % (strain))
				report_file.write(strain+' contained %i scaffolds in the REORDERED file, and %i in the original assembly\n' % (len(scaffolds_reorderedFile),len(scaffolds_originalAssembly)))
				report_file.write('====================Scaffolds not in reordered============================\n')
				scaffolds_notReordered =[]
				for scaffold in scaffolds_originalAssembly: #Check the scaffolds that are not in reordered
					if scaffold.name not in scaffolds_reorderedFile:
						scaffolds_notReordered.append(scaffold.name)
						report_file.write(scaffold.name+'\n')
				records = (r for r in SeqIO.parse(line.split('\t')[0], "fasta") if r.id in scaffolds_notReordered)
				count_notInReordered = SeqIO.write(records, current_folder+strain + '_notInReordered.fasta', "fasta")
		print "Scaffolds not reordered are saved in %s" % current_folder+strain + '_notInReordered.fasta'
		if count_notInReordered > 0:
			os.system(EMBOSS_info(current_folder+strain + '_notInReordered.fasta',strain+ '_notInReordered'))
		os.system(EMBOSS_info(current_folder+strain + '_reordered.fasta',strain+ '_reordered'))
		report_file.close()
		os.system(all_contigs(strain))

print "DONE!"
