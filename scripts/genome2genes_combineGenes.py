__author__ = 'Peris'

#!/usr/bin/env python
#coding: utf-8

import argparse
import os
import shutil
import glob
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import GC

helptext="""
This script is to extract genes and proteins from a genome assembly and aligned the orthologous genes based on table information
Authors: David Peris UW-Madison, Dept Genetics
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A table with StrainNames and presence of genes", type = str, default = None)
parser.add_argument("-f","--folder", help="Folder where the _indivGenes_renamed folders are stored", type = str, default = None)
parser.add_argument("-o","--outfolder", help="Folder where the combined sequences are stored", type = str, default = None)


args = parser.parse_args()

table_genes = open(args.input, 'r')

counter = 0

if not os.path.exists(args.folder+args.outfolder):
	os.makedirs(args.folder+args.outfolder)
	if not os.path.exists(args.folder+args.outfolder+'/Sequences_together'):
		os.makedirs(args.folder+args.outfolder+'/Sequences_together')

output_folder = args.folder+ args.outfolder + '/Sequences_together/'

for line in table_genes:
	if counter == 0:
		list_strains = line.strip().split('\t')[1:]
		counter += 1
	else:
		GeneName = line.split('\t')[0]
		LineValues = line.split('\t')[1:]
		LineValues = "\t".join(LineValues)
		if not '0' in LineValues:
			print GeneName
			catDNA_cmd = 'cat '
			catAA_cmd =  'cat '
			for StrainName in list_strains:
				catDNA_cmd += args.folder+StrainName+'_indivGenes_renamed/'+GeneName+'.dna.fasta '
				catAA_cmd += args.folder+StrainName+'_indivGenes_renamed/'+GeneName+'.prot.fasta '
			catDNA_cmd += '> ' + output_folder + GeneName+'.dna.fasta'
			catAA_cmd += '> ' + output_folder + GeneName+'.prot.fasta'
			#print catDNA_cmd
			#print catAA_cmd
			os.system(catDNA_cmd)
			os.system(catAA_cmd)
			counter += 1

print "%i Genes have been combined" % (counter)

print "DONE!"
