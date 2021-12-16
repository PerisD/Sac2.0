__author__ = 'Quinn'

import sys
import re #to do regular expressions
from numpy import * #to do the transposition
import numpy as np
import csv
import random


vcfName = sys.argv[1] #vcf file name
structureName = sys.argv[2] #just the prefix "all_SNPs.import" or "10k_SNPs.import" will be added to it later

strainList = list()
vcfArray = list()
simSNParray = list()
SNPinfo = list()
finalArray = list()

vcf = open(vcfName, 'r')
lines = vcf.readlines()
header = str(filter(lambda x:'#CHROM' in x, lines))
strainList = re.findall('tibet|yH\w{2,3}\d{1,3}', header)
#print(strainList)
strainList.insert(0, '\t')
print strainList

interval = int((len(lines)-1)/10000) #this will be the interval for even spacing of a subset of 10K SNPs for a simplified
start = random.randint(1,10000)
print(start)
steps = range(start, len(lines), interval)
for i in steps:
    currentLine = lines[i].strip('\n')
    SNPline = currentLine.split() #Split based on tabs
    counter = 9     #9th position is where genotype data starts
    if len(SNPline[3]) == 1 and len(SNPline[4]) == 1: #Check that they arn't indels
        SNPinfo = list()
        SNPid = SNPline[0]+":"+SNPline[1]  #concatnate chr and position for SNP name
        SNPinfo.append(SNPid)   #add it to the list
        while counter<len(SNPline):
            if re.match('\.', SNPline[counter]): #build the table with the SNP info
                SNPinfo.append(0)
                SNPinfo.append(0)
            elif re.match('0', SNPline[counter]):
                SNPinfo.append(0)
                SNPinfo.append(1)
            elif re.match('^1', SNPline[counter]):
                SNPinfo.append(1)
                SNPinfo.append(1)
            counter = counter + 1
    simSNParray.append(SNPinfo)


simSNPnparray = array(simSNParray)
simTransposed = simSNPnparray.T
simListTrans = simTransposed.tolist()

#could make two for loops going through each list
simListTrans[0].insert(0, strainList[0])
strainCounter = 1
snpCounter = 1
while strainCounter<len(strainList):
    simListTrans[snpCounter].insert(0, strainList[strainCounter])
    simListTrans[snpCounter+1].insert(0, strainList[strainCounter])
    strainCounter = strainCounter + 1
    snpCounter = snpCounter + 2

simOutStrName = structureName + "_10k_SNPs.str"


with open(simOutStrName, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ', lineterminator='\n')
    writer.writerows(simListTrans)
