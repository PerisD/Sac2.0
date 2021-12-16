#!/usr/bin/env python

import itertools
import pandas as pd
import glob
import argparse

helptext="""
This script will generate a csv table with a relation between StrainName and Haplotypes inferred from DnaSP v5
Authors: Peris IATA-CSIC, Paterna, Valencia, Spain
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--inputFile", help="A nexus haplotype file", type = str, default = None)
parser.add_argument("-o","--outputFile", help="A csv file Strain - Haplotype designation", type = str, default = None)

parser.set_defaults(spades=True)

args = parser.parse_args()

haplotypes_nex = open(args.inputFile,"r")
counter1 = 0
counter2 = 1
counter3 = 1
haplotypes_list = []
strains_list = []
for line in haplotypes_nex:
	if "DIMENSIONS NTAX" in line:
		number_haplotypes = line.strip().split('=')[-1][:-1]
		number_haplotypes = int(number_haplotypes)
		print number_haplotypes
	if "[Hap#  Freq. Sequences]" in line:
		counter1 += 1
	if counter1 == 2:
		if counter2 == 1:
				counter2 += 1
		elif counter2 >= 2 and counter3 <= number_haplotypes:
				line = line.strip().split(':')
				haplotype_now = line[0]
				haplotype_now = haplotype_now[1:]
				strains = line[1].split()
				number_strains_now = strains[0]
				strains = strains[1:]
				strains[-1] = strains[-1][:-1]
				counter3 += 1
				print haplotype_now
				print number_strains_now
				print strains
				total_haplotypes_now = list(itertools.repeat(haplotype_now,int(number_strains_now)))
				counter4 = 0
				while counter4 < len(total_haplotypes_now):
					haplotypes_list.append(total_haplotypes_now[counter4])
					counter4 += 1
				for strain_name in strains:
					strains_list.append(strain_name)

haplotypes_strains = pd.DataFrame({'StrainName': strains_list,'Haplotype': haplotypes_list})
haplotypes_strains.to_csv(args.outputFile+'.csv', index=False, sep=',')