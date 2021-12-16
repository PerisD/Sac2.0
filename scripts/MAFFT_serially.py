#! /usr/bin/env python
#coding: utf-8

import os
import glob
import argparse

helptext="""
This script will rename the sequence names if in BUSCO was added the run_ name
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="PATH to the folder where alignments are stored", type = str, default = None)
parser.add_argument("-o","--output", help="Output folder of the new fasta files", type = str, default = "fasta_aligned")
parser.add_argument("-t","--threads", help="number of threads", type = str, default = "4")

args = parser.parse_args()

input_folder = args.input
output_folder = args.output

def mafft_serially(fasta,new_fasta):
	cmdA = 'mafft --op 1.53 --ep 0.123 --auto --thread '
	cmdA += args.threads + ' ' 
	cmdA += fasta + ' > ' + output_folder + '/' + new_fasta
	return cmdA

if os.path.exists(output_folder):
    print "Folder %s is not necessary" % (output_folder)
else:
    os.makedirs(output_folder)

list_fastas = glob.glob(input_folder+'*.fas')
total_files = len(list_fastas)

count = 0
for fasta in list_fastas:
	new_fasta = fasta.split('/')[-1]
	print "Aligning %s: %i/%i" % (new_fasta,count+1,total_files)
	os.system(mafft_serially(fasta,new_fasta))
	count += 1

print "Aligned of %i files is done" % count