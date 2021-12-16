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
parser.add_argument("-o","--output", help="Output folder of the new fasta files", type = str, default = "renamed_fastas")

args = parser.parse_args()

input_folder = args.input
output_folder = args.output

if os.path.exists(output_folder):
    print "Folder %s is not necessary" % (output_folder)
else:
    os.makedirs(output_folder)

list_fastas = glob.glob(input_folder+'*.fas')
total_files = len(list_fastas)

count = 0
for file in list_fastas:
	new_file = file.split('/')[-1]
	file_now = open(file, 'r')
	new_file_now = open(output_folder+'/'+new_file,'w')
	print "Renaming %s: %i/%i" % (new_file,count+1,total_files)
	for line in file_now:
		if line.startswith('>'):
			strain_name = line.split('_')[-1]
			line = ">" + strain_name
			new_file_now.write(line)
		else:
			new_file_now.write(line)
	new_file_now.close()
	file_now.close()
	count += 1

print "Rename of %i files is done" % count