__author__ = 'Quinn'
from Bio import SeqIO
import sys

strainListName = sys.argv[1]
outName = sys.argv[2]


outFileName = outName + ".fasta"

strainDict = {}
outFile = open(outFileName, 'w')
strainList = open(strainListName, 'r')
lines = strainList.readlines()
for line in lines :
	strain = line.strip()
	strainFasta = open(strain, 'r')
	strainString = ""
	for seq_record in SeqIO.parse(strainFasta, "fasta"):
			strainString = strainString + str(seq_record.seq)
	strainDict[strain] = strainString
	strainName = strain.split("/")[-1]
	strainName = strainName.split(".")[0]
	outFile.write(">" + strainName + "\n" + strainString + "\n")

