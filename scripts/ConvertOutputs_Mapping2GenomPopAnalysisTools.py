__author__ = 'Peris based on Quinn\'s scripts'

#!/usr/bin/env python
#coding: utf-8

import argparse
from subprocess import call
import os
import shutil
import glob
from Bio import SeqIO
import sys
import re #to do regular expressions
from numpy import * #to do the transposition
import numpy as np
import csv
import random

helptext="""
This script is to generate output files for being use in STRUCTURE
Authors: Based on Quinn Langdon scripts modified by Peris UW-Madison, Dept Genetics & IATA-CSIC
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A tab tabulated text file with StrainName and read paths used in Mapping4SNPs_v2.0.py, the reads will not be used BTW", type = str, default = None)
parser.add_argument("-t","--threads", help="Number of CPUs to be used", type = str, default = "8")
parser.add_argument("-r","--reference", help="Reference Genome path to map the reads", type = str, default = None)
parser.add_argument("-o","--output", help="the output name", type = str, default = None)
parser.add_argument("-s","--outStrain", help="the strain used as outgroup", type = str, default = None)
parser.add_argument("-f","--folder", help="folder where Mapping4SNPs_v2.0.py results are stored where strain folders are found", type = str, default = None)
parser.add_argument("-W","--windowSNPs", help="Window size to get the SNPs, default is 10000 bp", type = int, default = 10000)
parser.add_argument("-R","--RemoveStrain", help="strains to be removed from the input file, example to remove 2 strains: yHAB336", type = str, default = '')


parser.set_defaults(spades=True)

args = parser.parse_args()

info_path = open(args.input, 'r')

GenomePath = args.reference

#Step to combine variant.vcfs using GATK

list_of_strains = []
for line in info_path:
	line = line.split('\t')
	StrainName = line[0]
	if not StrainName == args.outStrain:
		list_of_strains.append(StrainName)

#print "Strain list before removing:"+'_'.join(list_of_strains)

strains2remove = args.RemoveStrain
output_tagRm = ''

if strains2remove != '':
	strains2remove = strains2remove.split(',')
	counterRmStrains = 0
	for strain in strains2remove:
		print "Strain being removed:"+strain
		list_of_strains.remove(strain)
		counterRmStrains += 1
	output_tagRm = '_' + str(counterRmStrains) + 'filtered'

print list_of_strains

info4Structure = open(args.output+'_Info4Structure.txt','w')
info4Structure.write('Number of individuals:'+str(len(list_of_strains))+'\n')

gatk_cm1 = '/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/gatk3/GenomeAnalysisTK.jar -T CombineVariants -R ' + GenomePath + ' '

for Strain in list_of_strains:
	print "%s calling variant" % (Strain)
	gatk_cm1 += '--variant ' + args.folder + Strain + '/' + Strain + '_SNP/' + Strain + '_variants.vcf '

gatk_cm1 += '-o ' + args.output + '_merged' + output_tagRm +'.vcf -genotypeMergeOptions UNIQUIFY'
print "GATK:"+gatk_cm1
os.system(gatk_cm1)

vcfArray = list()
simSNParray = list()
SNPinfo = list()
finalArray = list()

vcf = open(args.output + '_merged' + output_tagRm +'.vcf', 'r')
lines = vcf.readlines()
for VCFline in lines:
	if "#CHROM" in VCFline:
		VCFline_strains = VCFline.split('\t')[9:]
#print VCFline
list_strains_vcf = []
for i in VCFline_strains:
	i = i.split('.')[0]
	list_strains_vcf.append(i)

unusedHeader = filter(lambda x:'##' in x, lines)
main = filter(lambda x:re.match('chrI|chrV|chrX',x), lines)
header = str(filter(lambda x:'#CHROM' in x, lines))

list_strains_vcf.insert(0,'\t')
print list_strains_vcf

interval = int((len(lines)-1)/args.windowSNPs) #this will be the interval for even spacing of a subset of 10K SNPs for a simplified
start = random.randint(1,interval)
print(start)
steps = range(start, len(main), interval)
print("There should be "+str(len(steps))+" SNPs per line")

for i in steps:
    currentLine = main[i].strip('\n')
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

info4Structure.write('Number of SNPs pulled out:'+str(len(steps)))
info4Structure.close()

simSNPnparray = array(simSNParray)
simTransposed = simSNPnparray.T
simListTrans = simTransposed.tolist()

#could make two for loops going through each list
simListTrans[0].insert(0, list_strains_vcf[0])
strainCounter = 1
snpCounter = 1
while strainCounter<len(list_strains_vcf):
	simListTrans[snpCounter].insert(0, list_strains_vcf[strainCounter])
	simListTrans[snpCounter+1].insert(0, list_strains_vcf[strainCounter])
	strainCounter = strainCounter + 1
	snpCounter = snpCounter + 2

simOutStrName = args.output + output_tagRm + "_"+ str(args.windowSNPs/1000)+"K_SNPs.str"

with open(simOutStrName, 'w') as csvfile:
	writer = csv.writer(csvfile, delimiter=' ', lineterminator='\n')
	writer.writerows(simListTrans)

print "Done!"
