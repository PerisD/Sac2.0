#! /usr/bin/env python
#coding: utf-8

import os
import argparse
import glob
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

helptext="""
This script will check the length of the BUSCO sequences for a specific specimen
Authors: David Peris IATA-CSIC
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--inputFolder", help="folder where trimmed fasta are stored", type = str, default = None)
parser.add_argument("-o","--originalFolder", help="folder where original fasta are stored", type = str, default = None)
parser.add_argument("-a","--alnextension", help="how your original sequences are ending, default _nt_aln.fas", type = str, default = "_nt_aln.fas")
parser.add_argument("-t","--trextension", help="how your trimmed sequences are ending, default _nt_tr.fas", type = str, default = "_nt_tr.fas")


parser.set_defaults()

args = parser.parse_args()

extension1 = args.alnextension
fasta_files = glob.glob(args.originalFolder+"*" + extension1)
extension2 = args.trextension

info_fasta = open("trimmed_lengthInfo.txt",'w')
info_fasta.write('\t'.join(["GeneName","Original Length","Trimmed Lenght","Percentage trimmed\n"]))
list_notrimmed = open("genes_absentfromtrimmedFolder.txt",'w')
counter1 = 0
for iFasta in fasta_files:
	counter1 += 1
	FastaName = iFasta.split('/')[-1].split(extension1)[0]
	#print FastapName + ": " + str(counter1) + " of " + str(len(fasta_files))
	for record1 in SeqIO.parse(iFasta,'fasta'):
		length_fastaAln = len(str(record1.seq).replace("-",""))
		if os.path.isfile(args.inputFolder+FastaName+extension2):
			for record2 in SeqIO.parse(args.inputFolder+FastaName+extension2,'fasta'):
				length_fastaTr = len(str(record2.seq).replace("-",""))
				break
			info_fasta.write('\t'.join([FastaName,str(length_fastaAln),str(length_fastaTr),str(((float(length_fastaAln)-float(length_fastaTr))/float(length_fastaAln))*100)])+"\n")
		else:
			list_notrimmed.write(FastaName+"\n")
		break
info_fasta.close()
list_notrimmed.close()

length_df = pd.read_csv("trimmed_lengthInfo.txt", sep = "\t")
_ = plt.hist(length_df['Percentage trimmed'])
_ = plt.xlabel('percent of trimmed length')
_ = plt.ylabel('number of genes')
plt.savefig("distribution_trimmed_genes.pdf")

print "A total of " + str(counter1) + " sequences have been checked!\nDone!"