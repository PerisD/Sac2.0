__author__ = 'Peris'

#! /usr/bin/env python
#coding: utf-8

import sys,os
import argparse

helptext="""
Script for trimming using trimmomatic automatically
Authors: David Peris UW-Madison, Dept Genetics & IATA-CSIC
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="folder where read files to be trimmed are located", type = str, default = None)
parser.add_argument("-o","--output", help="Folder to save the trimmed files", type = str, default = None)

parser.set_defaults(spades=True)

args = parser.parse_args()

input_folder = args.input
output_folder = args.output + "/"

if not os.path.exists(output_folder):
	os.makedirs(output_folder)

#We will create a list of our files in the input folder
list_of_files = 0
list_of_files = os.listdir(input_folder)

def trimmomatic(fasta_file1,fasta_file2): 
	cmdA = 'java -jar /home/GLBRCORG/pnavarro/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 2 -phred33 -trimlog analytics.txt ' + input_folder
	cmdA += fasta_file1
	cmdA += ' ' + input_folder
	cmdA += fasta_file2
	cmdA += ' ' +output_folder+fasta_file1.split('.')[0]+'.fq ' +output_folder+ fasta_file1+'.tm-se.fq ' +output_folder+ fasta_file2.split('.')[0]+'.fq ' +output_folder+ fasta_file2+'.tm-se.fq '
	cmdA += 'ILLUMINACLIP:/home/GLBRCORG/pnavarro/software/Trimmomatic-0.33/adapters/TruSeq2-PE.fa:2:30:10 TRAILING:3 MINLEN:25'
	print cmdA
	return cmdA

list_of_files = sorted(list_of_files)
length = len(list_of_files) / 2
counter = 1
while counter <= length:
	read_name = list_of_files[0].split('.')[0]
	strain_name = read_name.split('_')[0]
	if os.path.isfile(input_folder+strain_name+'_1.fastq'):
		fasta_file1 = strain_name+'_1.fastq'
		fasta_file2 = strain_name+'_2.fastq'
		reads_to_remove = [fasta_file1,fasta_file2]
		for objects in reads_to_remove:
			list_of_files = filter(lambda a: a !=objects, list_of_files)
		counter += 1
		os.system(trimmomatic(fasta_file1,fasta_file2))
	else:
		fasta_file1 = strain_name+'_1.fq'
		fasta_file2 = strain_name+'_2.fq'
		reads_to_remove = [fasta_file1,fasta_file2]
		for objects in reads_to_remove:
			list_of_files = filter(lambda a: a !=objects, list_of_files)
		counter += 1
		os.system(trimmomatic(fasta_file1,fasta_file2))

print "done"
