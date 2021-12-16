#! /usr/bin/env python
#coding: utf-8

import sys,os
import argparse

helptext="""
Python 2 script
This script will make nuc blastdb automatically.
Authors: David Peris University of Oslo
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
#parser.add_argument("-i","--input", help="PATH to the folder where assemblies to be used are stored", type = str, default = None)
parser.add_argument("-a","--assemblies", help="assembly fasta file or list of assemblies in a text file .txt", type = str, default = None)
parser.add_argument("-o","--output", help="Folder to store the final blast databases", type = str, default = None)

parser.set_defaults()

args = parser.parse_args()

#input_folder = args.input
output_folder = args.output

#We will create a list of our files in the current path
list_of_files = []
if "txt" in args.assemblies.split(".")[-1]:
	assembly_list = open(args.assemblies,'r')
	for iLine in assembly_list:
		list_of_files.append(iLine.strip())
else:
	list_of_files.append(args.assemblies)

for fasta_file in list_of_files: # runs the makeblast command for all genomes 
	strain_name = fasta_file.split('/')[-1].split('.')[0]
	cmdA = 'makeblastdb -in ' + fasta_file
	cmdA += ' -dbtype nucl -out ' + output_folder + "/"
	cmdA += strain_name
	print cmdA
	os.system(cmdA)

print "done"
