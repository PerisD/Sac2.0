__author__ = 'Quinn modified by Peris'

from Bio import SeqIO
import sys
import re
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
import time
import pandas as pd
import argparse

helptext="""
This script is to generate output files for being use in PCAdmix
Authors: Based on Quinn Langdon scripts modified by Peris UW-Madison, Dept Genetics & IATA-CSIC
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="fasta file generated from Mapping4SNPs", type = str, default = None)
parser.add_argument("-p","--popName", help="Name of the list of strains in particular pop to generate the next files, i.e. Admix for a file strainList_Admix.txt", type = str, default = None)
parser.add_argument("-c","--chromosomeInfo", help="A tabulated file generated during maskCov_FASTA_byStrain-d_v2.1, pairwiseDivFile_chromosome.txt, with size of each chromosome and accumulated size", type = str, default = None)
parser.add_argument("-s","--strainInfo", help="Strain names and population information", type = str, default = None)
parser.add_argument("-r","--recombination", help="Recombination rate for each chromosome, we are using the Scer this might change in the future", type = str, default = None)


parser.set_defaults(spades=True)

args = parser.parse_args()

fastaName = args.input
ignoreN = float(0.0)
popName = args.popName
chromosomeInfo = args.chromosomeInfo
strainInfo_csv = args.strainInfo
recombination_csv =  args.recombination

popFileName = "strainList_"+popName+".txt"
popFile = open(popFileName, 'r')
popStrains = popFile.read().splitlines()
if re.match('chr', fastaName):
    fastaNameSplit = fastaName.split('.')
    chrName = fastaNameSplit[0]
    popOutFileName = popName+"_"+chrName+".txt"
else:
    popOutFileName = popName+".txt"
popOutFile = open(popOutFileName, 'w')
popOutFile.write("I\tsnpID")

for strain in popStrains:
    popOutFile.write("\t"+strain+"\t"+strain+"\t")

#Section to get the information about the chromosome and length
chromosomeInfo_df = pd.read_csv(chromosomeInfo, sep='\t')
chrList = chromosomeInfo_df['Chr'].tolist()

counterA = 0
maxSizeChr_list = []
sizeChr_accumulated = 0
for i in chrList:
    if counterA == 0:
        sizeChr_temp =  chromosomeInfo_df.loc[chromosomeInfo_df['Chr'] == i,'Size'].iloc[0]
        globals()[i+"max"] = sizeChr_temp
        sizeChr_accumulated += sizeChr_temp
        maxSizeChr_list.append(i+"max")
    if counterA >= 1:
        sizeChr_temp =  chromosomeInfo_df.loc[chromosomeInfo_df['Chr'] == i,'Size'].iloc[0]
        sizeChr_accumulated += sizeChr_temp
        globals()[i+"max"] = sizeChr_accumulated
        maxSizeChr_list.append(i+"max")
    counterA +=1

chrPosDict = {}
byChrPos = list()
for chr in chrList:
    chrPosDict[chr] = list()


strainDict = {}
strainList = list()
length = int()
fasta = open(fastaName, 'r')
for seq_record in SeqIO.parse(fasta, "fasta"):
    strainList.append(seq_record.id)
    strainDict[seq_record.id] = seq_record.seq
    length = len(seq_record.seq)

recombination_csv_f = open(recombination_csv,'r')
recombination_csv_f.next()
recomboDict = {}
for line in recombination_csv_f:
    line = line.strip().split(',')
    recomboDict[line[0]] = line[1]

strainInfo_csv_f = open(strainInfo_csv,'r')
strainInfo_csv_f.next()
popDict = {}
for line in strainInfo_csv_f:
    line = line.strip().split(',')
    popDict[line[0]] = line[2].replace(' ','_')

snpDict = {}
for strain in strainList:
    snpDict[strain] = str()

SNPpos = str()

first = strainList[0]
refMissing = 0
hasData = 0

start = time.time()
numSNPs = 0
for i in xrange(0, length): #go through the genome position by position
    compare = strainDict[first][i]
    if re.match('A|T|C|G', compare):
        hasData += 1
        matches = 0
        missing = 0
        SNP = 0
        propMiss = 0.0
        if re.match('N', compare):
            missing = 1
        for strain in strainList[1: len(strainList)]:
            if re.match(compare, strainDict[strain][i]):
                matches += 1
            else:
                if re.match('A|T|C|G', strainDict[strain][i]):
                    SNP += 1
                else:
                    missing += 1
        propMiss = float(missing)/len(strainList)
        if SNP>0 and propMiss<=ignoreN:
            numSNPs += 1
            popOutFile.write("\n")
            endTime = time.time()-start
            if 60 < endTime < 3600:
                min = int(endTime)/60
                sec = int(endTime-(min*60))
                elapsedTime = str(min) + " mins " + str(sec) + " secs"
            elif 3600 < endTime < 86400:
                hr = int(endTime)/3600
                min = int((endTime - (hr*3600))/60)
                sec = int(endTime - ((hr*60)*60 + (min*60)))
                elapsedTime = str(hr) + " hrs " + str(min) + " mins " + str(sec) + " secs"
            elif 86400 < endTime < 604800:
                day = int(endTime)/86400
                hr = int((endTime-(day*86400))/3600)
                min = int((endTime - (hr*3600+day*86400))/60)
                sec = int(endTime - ((day*86400) + (hr*3600) + (min*60)))
                elapsedTime = str(day)  + " days " + str(hr) + " hrs " + str(min) + " mins " + str(sec) + " secs"
            else:
                elapsedTime = str(int(endTime)) + " secs"
            print("%.5f" % ((i/11633661.0)*100) + "%" + "\tElapsed time: " + elapsedTime)
            SNPid = str()
            truPos = i+1
            if re.match('chr', fastaName):
                SNPid = str(truPos)
                byChrPos.append(truPos)
            elif 0 < truPos < chrImax:
                chrPos = str(truPos)
                SNPid = "chrI:" + str(truPos)
                chrPosDict[chrList[0]].append(truPos)
            elif chrImax<truPos<chrIImax:
                chrPos = str(truPos-chrImax)
                SNPid = "chrII:" + str(truPos-chrImax)
                chrPosDict[chrList[1]].append(truPos)
            elif chrIImax < truPos < chrIIImax:
                chrPos = str(truPos-chrIImax)
                SNPid = "chrIII:" + str(truPos-chrIImax)
                chrPosDict[chrList[2]].append(truPos)
            elif chrIIImax < truPos < chrIVmax:
                chrPos = str(truPos-chrIIImax)
                SNPid = "chrIV:" + str(truPos-chrIIImax)
                chrPosDict[chrList[3]].append(truPos)
            elif chrIVmax < truPos < chrVmax:
                chrPos = str(truPos-chrIVmax)
                SNPid = "chrV:" + str(truPos-chrIVmax)
                chrPosDict[chrList[4]].append(truPos)
            elif chrVmax < truPos < chrVImax:
                chrPos = str(truPos-chrVmax)
                SNPid = "chrVI:" + str(truPos-chrVmax)
                chrPosDict[chrList[5]].append(truPos)
            elif chrVImax < truPos < chrVIImax:
                chrPos = str(truPos-chrVImax)
                SNPid = "chrVII:" + str(truPos-chrVImax)
                chrPosDict[chrList[6]].append(truPos)
            elif chrVIImax < truPos < chrVIIImax:
                chrPos = str(truPos-chrVIImax)
                SNPid = "chrVIII:" + str(truPos-chrVIImax)
                chrPosDict[chrList[7]].append(truPos)
            elif chrVIIImax < truPos < chrIXmax:
                chrPos = str(truPos-chrVIIImax)
                SNPid = "chrIX:" + str(truPos-chrVIIImax)
                chrPosDict[chrList[8]].append(truPos)
            elif chrIXmax < truPos < chrXmax:
                chrPos = str(truPos-chrIXmax)
                SNPid = "chrX:" + str(truPos-chrIXmax)
                chrPosDict[chrList[9]].append(truPos)
            elif chrXmax < truPos < chrXImax:
                chrPos = str(truPos-chrXmax)
                SNPid = "chrXI:" + str(truPos-chrXmax)
                chrPosDict[chrList[10]].append(truPos)
            elif chrXImax < truPos < chrXIImax:
                chrPos = str(truPos-chrXImax)
                SNPid = "chrXII:" + str(truPos-chrXImax)
                chrPosDict[chrList[11]].append(truPos)
            elif chrXIImax < truPos < chrXIIImax:
                chrPos = str(truPos-chrXIImax)
                SNPid = "chrXIII:" + str(truPos-chrXIImax)
                chrPosDict[chrList[12]].append(truPos)
            elif chrXIIImax < truPos < chrXIVmax:
                chrPos = str(truPos-chrXIIImax)
                SNPid = "chrXIV:" + str(truPos-chrXIIImax)
                chrPosDict[chrList[13]].append(truPos)
            elif chrXIVmax < truPos < chrXVmax:
                chrPos = str(truPos-chrXIVmax)
                SNPid = "chrXV:" + str(truPos-chrXIVmax)
                chrPosDict[chrList[14]].append(truPos)
            elif chrXVmax < truPos < chrXVImax:
                chrPos = str(truPos-chrXVmax)
                SNPid = "chrXVI:" + str(truPos-chrXVmax)
                chrPosDict[chrList[15]].append(truPos)
            elif chrXVImax < truPos < chr2_micronmax:
                chrPos = str(truPos-chrXVImax)
                SNPid = "chr2_micron:" + str(truPos-chrXVImax)
                chrPosDict[chrList[16]].append(truPos)
            elif chr2_micronmax < truPos < chrMTmax:
                chrPos = str(truPos-chr2_micronmax)
                SNPid = "chrMT:"+chrPos
                chrPosDict[chrList[17]].append(truPos)
            elif truPos > chrMTmax:
                chr = "unknown"
                chrPos = str(truPos-chrMTmax)
                SNPid = "unknown:" + str(truPos-chrMTmax)
            #SNPid = chr+":"+chrPos
            popOutFile.write("M\t"+SNPid)
            SNPpos =  SNPpos + " " + str(truPos)
            for strain in popStrains:
                strainSNP = strainDict[strain][i]
                popOutFile.write("\t"+strainSNP+"\t"+strainSNP)
    else:
        refMissing += 1
            #i += 1

