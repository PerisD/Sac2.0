__author__ = 'Quinn Langdon modified by Peris'

#!/usr/bin/env python
from Bio import SeqIO
import sys

#Python script to take list of fastas and create one file for each chromosome

strainListName = sys.argv[1]
chrListName = sys.argv[2]

chrList = open(chrListName, 'r')
lines = chrList.readlines()
for line in lines:
	chr = line.strip()
	outputFileName = chr + '.fasta'
	outFile = open(outputFileName, 'w')
	strainList = open(strainListName, 'r')
	lines = strainList.readlines()
	chrOutput = []
	for line in lines :
		strain = line.strip()
		strainFasta = open(strain, 'r')
		#to read in the fasta
		strainDict = {}
		for seq_record in SeqIO.parse(strainFasta, "fasta"):
				strainDict[seq_record.id] = seq_record.seq
		chrSequence = str(strainDict[chr])
		#strainFastaOutput = strain + "\n" + chrSequence
		#chrOutput.append(strain)
		#chrOutput.append(chrSequence)
		strainName = strain.split("/")[-1]
		strainName = strainName.split(".")[0]
		outFile.write(">" + strainName + "\n")
		outFile.write(chrSequence +"\n")