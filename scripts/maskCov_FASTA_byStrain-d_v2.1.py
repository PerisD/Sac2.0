__author__ = 'Quinn modfied by Peris'

#!/usr/bin/env python
from Bio import SeqIO
import sys
import re
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

#Python script to extract data from fast and vcf
#print 'Arguments: coverage quant by strain  \n'

strainCovQuantName = sys.argv[1]
chrListName = sys.argv[2]
fastaToChangeName = sys.argv[3]
refGenome = sys.argv[4]
pairwiseFile = 'pairwiseDivFile_chromosome.txt'
#Dictionary of chromosome lengths
#chrLen = {'chrI': 175444, 'chrII': 1274804, 'chrIII': 305615, 'chrIV': 982538, 'chrV': 588913, 'chrVI': 264757, 'chrVII': 1051730, 'chrVIII': 741893, 'chrIX': 401362, 'chrX': 747934, 'chrXI': 632881, 'chrXII': 1034152, 'chrXIII': 953685, 'chrXIV': 768019, 'chrXV': 813827, 'chrXVI': 896107}
#comboChrLen = {'Seub_chrI': 175444, 'Seub_chrII': 1274804, 'Seub_chrIII': 305615, 'Seub_chrIV': 982538, 'Seub_chrV': 588913, 'Seub_chrVI': 264757, 'Seub_chrVII': 1051730, 'Seub_chrVIII': 741893, 'Seub_chrIX': 401362, 'Seub_chrX': 747934, 'Seub_chrXI': 632881, 'Seub_chrXII': 1034152, 'Seub_chrXIII': 953685, 'Seub_chrXIV': 768019, 'Seub_chrXV': 813827, 'Seub_chrXVI': 896107, 'Scer_chrI': 230218, 'Scer_chrII': 813184, 'Scer_chrIII': 316620, 'Scer_chrIV': 1531933, 'Scer_chrV': 576874, 'Scer_chrVI': 270161, 'Scer_chrVII': 1090940, 'Scer_chrVIII': 562643, 'Scer_chrIX': 439888, 'Scer_chrX': 745751, 'Scer_chrXI': 666816, 'Scer_chrXII': 1078177, 'Scer_chrXIII': 924431, 'Scer_chrXIV': 784333, 'Scer_chrXV': 1091291, 'Scer_chrXVI': 948066}

temp_length = []
pairwiseFile = open(pairwiseFile,'w')
pairwiseFile.write('Chr\tSize\tPosition_Genome_Start\n')

position = 0
for index, record in enumerate(SeqIO.parse(refGenome,"fasta")):
	chr_temp = record.id
	chr_size_temp = len(record.seq)
	line_now = [chr_temp,str(chr_size_temp),str(position)+'\n']
	position += chr_size_temp
	pairwiseFile.write('\t'.join(line_now))
	temp_length.append(len(record.seq))
genomeLen = float(sum(temp_length))
pairwiseFile.close()

maskingInfo = ""
strainCovQuant = open(strainCovQuantName, 'r')
strains = strainCovQuant.readlines()
strain = strains[1]
lowerCov = int(10)
strainLine = strain.strip('\n')
covInfo = strainLine.split()
strainName = re.sub(r'"','',covInfo[0])
lowerCheck = int(covInfo[1])
upperCov = int(float(covInfo[2]))
if lowerCheck < lowerCov:
    #print('Lower < 10')
    lowerCov = lowerCheck
maskingInfo = maskingInfo + strainName + ":\tlowLimit=" + str(lowerCov) + "\tupperLimit=" + str(upperCov) + "\n"
strainDict = {}
genomeDict = {}
#fastaToChangeName = strainName+".fasta"
fastaToChange = open(fastaToChangeName, 'r')
for seq_record in SeqIO.parse(fastaToChange, "fasta"):
    strainDict[seq_record.id] = seq_record.seq
    idStr = str(seq_record.id)
    seqStr = str(seq_record.seq)
    genomeDict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())
Ncount = 0
for key in genomeDict:
    Ncount = Ncount + genomeDict[key].count("N")
lenToMaskUpper=0
lenToMaskLower=0
existingN = 0
mito = 0
strainBed = strainName +".bedgraph"
bed = open(strainBed, 'r')
lines = bed.readlines()
for line in lines:
    currentLine = line.strip('\n')
    info = currentLine.split()
    chrom = info[0]
    chromstart = int(info[1])-1
    value = int(info[2])
    if value> upperCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
            genomeLen = 23704987.0
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskUpper = lenToMaskUpper + 1
    elif value<=lowerCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskLower = lenToMaskLower + 1
totalCovMasked = lenToMaskLower + lenToMaskUpper
totalMasked = totalCovMasked + Ncount - existingN
rawNumMasked =  "\t\tExisting N=" + str(existingN) + "\tLower to Mask=" + str(lenToMaskLower) + "\tUpper to Mask=" + str(lenToMaskUpper)

percentMasked =  "\t\tExisting N=" + "%.2f" % ((float(existingN)/genomeLen)*100) + "%\tLower to Mask=" + "%.2f" % ((float(lenToMaskLower)/genomeLen)*100) + "%\tUpper to Mask=" + "%.2f" % ((float(lenToMaskUpper)/genomeLen)*100) + "%"
strTotalMasked = "\t\tRaw Total Masked=" + str(totalMasked) + "\tPercent Total Masked=" + "%.2f" % ((totalMasked/genomeLen)*100) + "%"
#print(rawNumMasked)
#print(percentMasked)
#print(strTotalMasked)
maskingInfo = maskingInfo + rawNumMasked + "\n" + percentMasked + "\n" + strTotalMasked + "\n"

infoFileName = strainName+"_covMaskedInfo.txt"
infoFile = open(infoFileName, 'w')
infoFile.write(maskingInfo)
infoFile.close()

outputFileName = strainName+"_covMasked.fasta"
outFile = open(outputFileName, 'w')
#if re.match('.+lager', strainName):
#    chrListName = "chromosome_combo.txt"
#else:
#    chrListName = "chromosome.txt"
chrList = open(chrListName, 'r')
chrLines = chrList.readlines()
for chrs in chrLines:
    chr = chrs.strip()
    outFile.write(">"+ chr + "\n")
    outFile.write(str(genomeDict[chr]) + "\n")
outFile.close()