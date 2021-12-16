#! /usr/bin/env python
#coding: utf-8

import os
import argparse
import glob

helptext="""
This script will combine all sequences for the specimens ran during busco, at this point it takes all specimens.
Authors: David Peris IATA-CSIC
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
#parser.add_argument("-i","--input_file", help="text file with StrainName and PATH to the assembly used to run prepareBUSCO", type = str, default = None)
parser.add_argument("-d","--dataBase", help="database used in BUSCO, i.e. agaricomycetes_odb10", type = str, default = None)
parser.add_argument("-l","--fasta2Parse", help="Genes_in_All.csv file generate by busco_fullTable.R to get the multifasta for genes completed in all specimens", type = str, default = None)
parser.add_argument("-o","--outputFolder", help="Folder to store the new multifastas", type = str, default = None)
parser.add_argument("-s","--strainList", help="strains to parse for renaming sequences in a txt file or the strain \
name as XXX or separated by commas XXX,YYY", type = str, default = None)


parser.set_defaults()

args = parser.parse_args()

if not os.path.exists(args.outputFolder + "/"):
	os.makedirs(args.outputFolder + "/")

fasta_files = open(args.fasta2Parse,'r')
next(fasta_files)
num_lines = 0
for iFasta in fasta_files:
	num_lines += 1
fasta_files.close()

#We need to rename the files first
if ".txt" in args.strainList:
	list_Strains = open(args.strainList,'r')
else:
	list_Strains = []
	if "," in args.strainList:
		list_Strains = args.strainList.split(",")
	else:
		list_Strains.append(args.strainList)
for iStrain in list_Strains:
	StrainName = iStrain.strip()
	list_files2rename = glob.glob(StrainName+"/run_"+args.dataBase+ "/busco_sequences/single_copy_busco_sequences/*.f*a")
	for iFile in list_files2rename:
		fasta_file = open(iFile,'r')
		for line in fasta_file:
			if line.startswith('>'):
				old_line = line.strip()
				new_line = ">"+StrainName
		fasta_file.close()
		modify_iFile = "sed -i 's|"+old_line+"|"+new_line+"|g' " + iFile
		#print(modify_iFile)
		os.system(modify_iFile)

fasta_files = open(args.fasta2Parse,'r')
next(fasta_files)
counter = 0
for iFasta in fasta_files:
	counter += 1
	iFasta = iFasta.split(',')[0]
	print(iFasta + ": " + str(counter) + " of " + str(num_lines))
	pullNT_cmd = "cat */run_" + args.dataBase+ "/busco_sequences/single_copy_busco_sequences/" + iFasta + "*.fna > " + args.outputFolder + "/" + iFasta + "_nt.fas"
	pullAA_cmd = "cat */run_" + args.dataBase+ "/busco_sequences/single_copy_busco_sequences/" + iFasta + "*.faa > " + args.outputFolder + "/" + iFasta + "_aa.fas"
	os.system(pullNT_cmd)
	os.system(pullAA_cmd)

print("A total of " + str(counter) + " sequences have been renamed and pulled\nDone!")