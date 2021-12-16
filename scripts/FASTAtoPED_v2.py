__author__ = 'Quinn modified by Peris D'

from Bio import SeqIO
import sys
import re
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
import time
import pandas as pd
import argparse

helptext="""
This script is to generate output files in PED format
ADMIXTOOLS requires some modifications of the PED format, such as swap the columns in .pedsnp, generate a PEDind file, and the recombination
map requires to be in M/bp, and this script generates the map in cM/kbp
Authors: Based on Quinn Langdon scripts modified by Peris UW-Madison, Dept Genetics & IATA-CSIC
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="fasta file generated from Mapping4SNPs", type = str, default = None)
parser.add_argument("-o","--output", help="The prefix for all output files", type = str, default = None)
parser.add_argument("-p","--percent", help="The proportion of strains that we will tolerate to have an N position, 0.0 tolerates 0 percent", type = float, default = 0.0)
parser.add_argument("-r","--recombination", help="Recombination rate for each chromosome, we are using the Scer this might change in the future", type = str, default = None)
parser.add_argument("-s","--strainInfo", help="Strain names and population information", type = str, default = None)
parser.add_argument("-c","--chromosomeInfo", help="A tabulated file generated during maskCov_FASTA_byStrain-d_v2.1, pairwiseDivFile_chromosome.txt, with size of each chromosome and accumulated size", type = str, default = None)


parser.set_defaults()

args = parser.parse_args()


fastaName = args.input
pedFileName = args.output + ".ped"
mapFileName = args.output + ".snp"
ignoreN = args.percent
recombination_csv =  args.recombination
strainInfo_csv = args.strainInfo
chromosomeInfo = args.chromosomeInfo


strainDict = {}
strainList = list()
length = int()
fasta = open(fastaName, 'r')
for seq_record in SeqIO.parse(fasta, "fasta"):
    strainList.append(seq_record.id)
    strainDict[seq_record.id] = seq_record.seq
    length = len(seq_record.seq)

strainInfo_csv_f = open(strainInfo_csv,'r')
strainInfo_csv_f.next()
popDict = {}
for line in strainInfo_csv_f:
    line = line.strip().split(',')
    popDict[line[0]] = line[2].replace(' ','_')

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

recombination_csv_f = open(recombination_csv,'r')
recombination_csv_f.next()
recomboDict = {}
for line in recombination_csv_f:
    line = line.strip().split(',')
    recomboDict[line[0]] = line[1]

snpDict = {}
for strain in strainList:
    snpDict[strain] = str()

SNPpos = str()
morgans = float()

first = strainList[0]
refMissing = 0
hasData = 0

start = time.time()


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
            #posOutFile.write("\n"+str(i+))
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
            chr = "1"
            SNPid = str()
            truPos = i+1
            if 0 < truPos < chrImax:
                chr = "1"
                chrPos = str(truPos)
                SNPid = "chrI:" + str(truPos)
                morgans = truPos*recomboDict["chrI"]
            elif chrImax<truPos<chrIImax:
                chr = "2"
                chrPos = str(truPos-chrImax)
                SNPid = "chrII:" + str(truPos-chrImax)
                morgans = (truPos-chrImax) * recomboDict["chrII"]
            elif chrIImax < truPos < chrIIImax:
                chr = "3"
                chrPos = str(truPos-chrIImax)
                SNPid = "chrIII:" + str(truPos-chrIImax)
                morgans = (truPos-chrIImax) * recomboDict["chrIII"]
            elif chrIIImax < truPos < chrIVmax:
                chr = "4"
                chrPos = str(truPos-chrIIImax)
                SNPid = "chrIV:" + str(truPos-chrIIImax)
                morgans = (truPos-chrIIImax) * recomboDict["chrIV"]
            elif chrIVmax < truPos < chrVmax:
                chr = "5"
                chrPos = str(truPos-chrIVmax)
                SNPid = "chrV:" + str(truPos-chrIVmax)
                morgans = (truPos-chrIVmax) * recomboDict["chrV"]
            elif chrVmax < truPos < chrVImax:
                chr = "6"
                chrPos = str(truPos-chrVmax)
                SNPid = "chrVI:" + str(truPos-chrVmax)
                morgans = (truPos-chrVmax) * recomboDict["chrVI"]
            elif chrVImax < truPos < chrVIImax:
                chr = "7"
                chrPos = str(truPos-chrVImax)
                SNPid = "chrVII:" + str(truPos-chrVImax)
                morgans = (truPos-chrVImax) * recomboDict["chrVII"]
            elif chrVIImax < truPos < chrVIIImax:
                chr = "8"
                chrPos = str(truPos-chrVIImax)
                SNPid = "chrVIII:" + str(truPos-chrVIImax)
                morgans = (truPos - chrVIImax) * recomboDict["chrVIII"]
            elif chrVIIImax < truPos < chrIXmax:
                chr = "9"
                chrPos = str(truPos-chrVIIImax)
                SNPid = "chrIX:" + str(truPos-chrVIIImax)
                morgans = (truPos - chrVIIImax) * recomboDict["chrIX"]
            elif chrIXmax < truPos < chrXmax:
                chr = "10"
                chrPos = str(truPos-chrIXmax)
                SNPid = "chrX:" + str(truPos-chrIXmax)
                morgans = (truPos - chrIXmax) * recomboDict["chrX"]
            elif chrXmax < truPos < chrXImax:
                chr = "11"
                chrPos = str(truPos-chrXmax)
                SNPid = "chrXI:" + str(truPos-chrXmax)
                morgans = (truPos - chrXmax) * recomboDict["chrXI"]
            elif chrXImax < truPos < chrXIImax:
                chr = "12"
                chrPos = str(truPos-chrXImax)
                SNPid = "chrXII:" + str(truPos-chrXImax)
                morgans = (truPos - chrXImax) * recomboDict["chrXII"]
            elif chrXIImax < truPos < chrXIIImax:
                chr  = "13"
                chrPos = str(truPos-chrXIImax)
                SNPid = "chrXIII:" + str(truPos-chrXIImax)
                morgans = (truPos - chrXIImax) * recomboDict["chrXIII"]
            elif chrXIIImax < truPos < chrXIVmax:
                chr = "14"
                chrPos = str(truPos-chrXIIImax)
                SNPid = "chrXIV:" + str(truPos-chrXIIImax)
                morgans = (truPos - chrXIIImax) * recomboDict["chrXIV"]
            elif chrXIVmax < truPos < chrXVmax:
                chr = "15"
                chrPos = str(truPos-chrXIVmax)
                SNPid = "chrXV:" + str(truPos-chrXIVmax)
                morgans = (truPos - chrXIVmax) * recomboDict["chrXV"]
            elif chrXVmax < truPos < chrXVImax:
                chr = "16"
                chrPos = str(truPos-chrXVmax)
                SNPid = "chrXVI:" + str(truPos-chrXVmax)
                morgans = (truPos - chrXVmax) * recomboDict["chrXVI"]
            elif chrXVImax < truPos < chr2_micronmax:
                chr = "17"
                chrPos = str(truPos-chrXVImax)
                SNPid = "chr2_micron:"+str(truPos-chrXVImax)
                morgans = (truPos - chrXVImax) * recomboDict["chr2_micron"]
            elif chr2_micronmax < truPos < chrMTmax:
                chr = "18"
                chrPos = str(truPos-chr2_micronmax)
                SNPid = "chrMT"+chrPos
                morgans = (truPos - chr2_micronmax) * recomboDict["chrMT"]
            elif truPos > chrMTmax:
                chr = "unknown"
                chrPos = str(truPos-chrXVImax)
                SNPid = "unknown" + str(truPos-chrXVImax)
                morgans = (truPos - chrXVImax) * recomboDict["unknown"]
            #SNPid = chr+":"+chrPos
            SNPpos =  SNPpos + chr + "\t" + SNPid+ "\t" + str(morgans) + "\t" + str(truPos) + "\n"
            for strain in strainList:
                snpDict[strain] = snpDict[strain] + strainDict[strain][i]+ " " + strainDict[strain][i] + " "
    else:
        refMissing += 1
            #i += 1

posOutFile = open(mapFileName, 'w')
posOutFile.write(SNPpos)
posOutFile.close()


outFile = open(pedFileName, 'w')
for strain in strainList:
    outFile.write(popDict[strain]+ "\t" + strain + "\t0\t0\t0\t0\t" +snpDict[strain] + "\n")

