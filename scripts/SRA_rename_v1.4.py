#!/usr/bin/env python

#########################################################################################################################
#Script to rename SRA accessions based on a tab file		                                                            #
#By David Peris {Genetics Lab, Biotechnology Center UW-Madison, WI}														#
#																														#
#usage: python SRA_rename.py ARG1 ARG2 ARG3 ARG4 ARG5																	#
#																														#
#ARG1 folder where SRA are located. Remember SRA files are in the format SRA_1.fastq/fq SRA_2.fastq/fq					#
#ARG2 extension name (i.e. fastq), useful to list all those files with *.fastq											#
#ARG3 file with StrainName	SRAaccession \t #StrainName will be use to rename SRA and create the folder					#
#ARG4 where to locate the new folders and fq files																		#
#																														#
#########################################################################################################################

import sys,os
import shutil
import glob

input_folder = sys.argv[1]
extension_name = sys.argv[2]
strain_SRA = sys.argv[3]
output_folder = sys.argv[4]

reference_SRA = open(strain_SRA,'r')

files_to_list = input_folder + '*_1.' + extension_name
list_of_files = glob.glob(files_to_list)
count = 0
SRA_notfound = []

for line in reference_SRA:
	SRA_name = line.split('\t')[1].strip()
	Strain_name = line.split('\t')[0]
	illumina_read1 = SRA_name + '_1.' + extension_name
	if any(illumina_read1 in s for s in list_of_files):
		illumina_read2 = illumina_read1.replace('_1','_2')
		if os.path.isfile(input_folder+illumina_read2):
			os.makedirs(output_folder+Strain_name)
			shutil.copy(input_folder+illumina_read1,output_folder+Strain_name+'/'+Strain_name+'_1.fq')
			shutil.copy(input_folder+illumina_read2,output_folder+Strain_name+'/'+Strain_name+'_2.fq')
		else:
			os.makedirs(output_folder+Strain_name)
			shutil.copy(input_folder+illumina_read1,output_folder+Strain_name+'/'+Strain_name+'_1.fq')
	else:
		count =+ 1
		SRA_notfound.append(SRA_name)

print "A total of %i SRA files were not found" % (count)

for file_notFound in SRA_notfound:
	print file_notFound

print "Done!"