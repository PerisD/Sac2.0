#! /usr/bin/env python
#coding: utf-8

helptext="""
This script is for BLASTing against multiple databases
Authors: David Peris UW-Madison, Dept Genetics & IATA-CSIC
"""

import sys,os
import argparse

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="can be a fasta file or a text file with the names of fasta to parse", type = str, default = None)
parser.add_argument("-o","--output", help="is the path to your folder where outputfiles will be saved", type = str, default = None)
parser.add_argument("-b","--blast", help="is the file with blast databases names", type = str, default = None)
parser.add_argument("-d","--database", help="is the folder where blast databases are stored", type = str, default = None)

args = parser.parse_args()

input_file = args.input
output_folder = args.output
blast_list = args.blast
blast_folder = args.database

if not os.path.exists(output_folder):
	os.makedirs(output_folder)
strainFasta_outputs = output_folder + 'fastas_individual/'
if not os.path.exists(strainFasta_outputs):
	os.makedirs(strainFasta_outputs)

if 'fas' in input_file.split('.')[1]:
	fasta_file = input_file
	databases = open(blast_list)
	for yeast in databases:
		run_all_db = 'blastn -db ' + blast_folder
		run_all_db += yeast.strip()
		run_all_db += ' -query '
		run_all_db += fasta_file # query fasta file
		run_all_db += ' -out '+ output_folder + "/" +yeast.strip()+'_'+fasta_file.split('/')[-1].split('.')[0]+'.txt'
		run_all_db += ' -outfmt \'6 qseqid sseqid sseq\' -evalue 0.01'
		print("run_all_db:"+run_all_db)
		os.system(run_all_db)
		make_fasta1 = 'awk \'BEGIN { OFS = "\\n" } { print ">'
		make_fasta1 += yeast.strip()
		make_fasta1 += '_"$1, $3 }\' <' + output_folder + "/" + yeast.strip() +'_'+ fasta_file.split('/')[-1].split('.')[0] + '.txt >>' + output_folder + fasta_file.split('/')[-1].split('.')[0] + '.fas'
		print("make_fasta_AllTogether:"+make_fasta1)
		os.system(make_fasta1)
		make_fasta2 = 'awk \'BEGIN { OFS = "\\n" } { print ">'
		make_fasta2 += yeast.strip()
		make_fasta2 += '_"$1, $3 }\' <' + output_folder + "/" + yeast.strip() +'_'+ fasta_file.split('/')[-1].split('.')[0] + '.txt >>' + strainFasta_outputs + yeast.strip() + '_AllFastas.fas'
		print("make_fasta_Individual:"+make_fasta2)
		os.system(make_fasta2)
else:
	fasta_list = open(input_file)
	for gene_file in fasta_list:
		fasta_file = gene_file.strip()
		databases = open(blast_list)
		for yeast in databases:
			run_all_db = 'blastn -db ' + blast_folder
			run_all_db += yeast.strip()
			run_all_db += ' -query '
			run_all_db += fasta_file # query fasta file
			run_all_db += ' -out '+ output_folder + "/" + yeast.strip()+'_'+fasta_file.split('/')[-1].split('.')[0]+'.txt'
			run_all_db += ' -outfmt \'6 qseqid sseqid sseq\' -evalue 0.01'
			print("run_all_db:"+run_all_db)
			os.system(run_all_db)
			make_fasta1 = 'awk \'BEGIN { OFS = "\\n" } { print ">'
			make_fasta1 += yeast.strip()
			make_fasta1 += '_"$1, $3 }\' <' + output_folder + "/" + yeast.strip() +'_'+ fasta_file.split('/')[-1].split('.')[0] + '.txt >>' + output_folder + fasta_file.split('/')[-1].split('.')[0] + '.fas'
			print("make_fasta:"+make_fasta1)
			os.system(make_fasta1)
			make_fasta2 = 'awk \'BEGIN { OFS = "\\n" } { print ">'
			make_fasta2 += yeast.strip()
			make_fasta2 += '_"$1, $3 }\' <' + output_folder + "/" + yeast.strip() +'_'+ fasta_file.split('/')[-1].split('.')[0] + '.txt >>' + strainFasta_outputs + yeast.strip() + "_AllFastas.fas"
			print("make_fasta_Individual:"+make_fasta2)
			os.system(make_fasta2)

remove_guion = "sed -i 's|-|N|g' " + strainFasta_outputs + "*.fas"
os.system(remove_guion)
print "DONE!"