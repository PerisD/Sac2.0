#! /usr/bin/env python
#coding: utf-8

######################################################################################
#Script for download and convert to fastq SRA datasets serially.                     #
#Authors: David Peris UW-Madison, Dept Genetics                                      #
#Usage: python download_SRA_serially.py INPUT OUTPUTFOLDER YES/NO                    #
#                                                                                    #
#INPUT a SRA accession number or a text file with a list of SRAs                     #
#OUTPUTFOLDER the folder where your fastq will be saved                              #
#YES or NO if your input is a list or just an accession number                       #
######################################################################################

import sys,os

SRA_files = sys.argv[1]
output_folder = sys.argv[2]
list_file = sys.argv[3]

downloaded_path = '~/ncbi/public/sra/'

if list_file == "NO":
	SRA_list = []
	SRA_list.append(SRA_files)
else:
	SRA_list = open(SRA_files)

def prefetch(SRA_file): #It is downloaded into the directory user/ncbi/public/sra/
	cmdA = 'prefetch -v ' + SRA_file
	return cmdA

def convert_fastq(SRA_file,output_folder):
	cmdB = 'fastq-dump --outdir ' + output_folder
	cmdB += ' --split-files ' + downloaded_path + SRA_file + '.sra'
	return cmdB

for SRA_file in SRA_list:
	SRA_file = SRA_file.strip()
	os.system(prefetch(SRA_file))
	os.system(convert_fastq(SRA_file,output_folder))

print "SRA files downloaded" 
