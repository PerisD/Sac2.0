__author__ = 'Quinn modified by Peris'

import sys
import re #to do regular expressions


vcfName = sys.argv[1] #vcf file name
outName = sys.argv[2]

vcfArray = list()
SNParray = list()
SNPinfo = list()
finalArray = list()
genotypeArray = list()
chrLen = dict()
chrList = list()
chrLenCumu = dict()
heterozygouteCounter = 0
homozygouteAltCounter = 0
homozygouteRefCounter = 0


outFile = open(outName+"_genotype.txt", 'w')
outFile.write("chr\tchrPos\tgenomePos\tGenotype\n")

vcf = open(vcfName, 'r')
lines = vcf.readlines()
for line in lines:
	currentLine = line.strip('\n')
	if re.match('##contig', currentLine): #look of the line with chr length info
		equalSplit = currentLine.split("=")
		commaSplit = equalSplit[2].split(",")
		chrName = commaSplit[0].strip("\s")
		chrLenStr = int(equalSplit[3].strip(">"))
		chrLen[chrName] = chrLenStr
		chrList.append(chrName)
		#SNParray.append(strainList)
	elif re.match('^#CHROM', currentLine):
		counter=0
		culumative = 0
		#print(chrList)
		for chr in chrList:
			chrLenCumu[chr] = culumative
			culumative = culumative+int(chrLen[chr])
#			print culumative
	elif re.match('^\w', currentLine): #Look for the data lines, needs to be changed for other line starters
		SNPline = currentLine.split() #Split based on tabs
		counter = 9	 #9th position is where genotype data starts
		if len(SNPline[3]) == 1 and len(SNPline[4]) == 1: #Check that they aren't indels
			SNPinfo = list()
			SNPid = SNPline[0]+":"+SNPline[1]  #concatnate chr and position for SNP name
			chr = SNPline[0]
			chrPos = SNPline[1]
			outFile.write(chr+"\t")
			outFile.write(chrPos+"\t")
			genomePos = chrLenCumu[chr]+int(SNPline[1])
			outFile.write(str(genomePos)+"\t")
			SNPinfo.append(SNPid)   #add it to the list
			if re.match('\.', SNPline[9]): #build the table with the SNP info
				SNPinfo.append(1) #if homozygous for the ref
				outFile.write("1" + "\t")
				#print(counter - 9)
				#print(strainList[counter-9])
				homozygouteRefCounter = homozygouteRefCounter+1
			elif re.match('0', SNPline[9]):
				SNPinfo.append(2) #if heterozygous
				outFile.write("2" + "\t")
				heterozygouteCounter = heterozygouteCounter+1
			elif re.match('^1', SNPline[9]):
				SNPinfo.append(3) #if homozygous for the alternate allele
				outFile.write("3" + "\t")
				homozygouteAltCounter = homozygouteAltCounter+1
		#outFile.write(SNPinfo + "\n")
			outFile.write("\n")
		SNParray.append(SNPinfo)


genotypeArray = (homozygouteRefCounter, heterozygouteCounter, homozygouteAltCounter)

genotypeOutName = outName + "_count.txt"
genotypeOut = open(genotypeOutName, 'w')
genotypeOut.write("Homozygous Ref\tHeterozygous\tHomozygous Alt\tTotal SNPs\tPercent Het of SNPs\tPercent Het of Genome\n")
genotypeOut.write(str(homozygouteRefCounter)+"\t")
genotypeOut.write(str(heterozygouteCounter)+"\t")
genotypeOut.write(str(homozygouteAltCounter)+"\t")
total = homozygouteAltCounter+homozygouteRefCounter+heterozygouteCounter
genotypeOut.write(str(total)+"\t")
genotypeOut.write(str((float(heterozygouteCounter)/total)*100)+"\t")
genotypeOut.write(str((float(heterozygouteCounter)/float(culumative))*100)+"%\n")

