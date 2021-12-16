__author__ = 'Miguel Morard modified by Peris'

#!/usr/bin/env python
#coding: utf-8

import sys
import csv
import argparse
import os
import glob

helptext="""
This script is to convert VCF file to tab file and get the heterozigosity distribution of allele frequencies
Authors: Miguel Morard and David Peris IATA-CSIC
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
#I make the input as the folder where other folders from the mapping4SNPs script was run
parser.add_argument("-i","--input", help="A folder with mapping4SNPs was run", type = str, default = None)
#parser.add_argument("-o","--output", help="The new table input of the snp_distrib script", type = str, default = None)
parser.add_argument("-p","--chromos", help="Information about the length of the chromosomes", type = str, default = "pairwiseDivFile_chromosome.txt")
parser.add_argument("-f","--format", help="version of the vcf, choose or 4.1 (my mapping) or 4.2(Miguel mapping)", type = str, default = "4.1")


parser.set_defaults(format=True)

args = parser.parse_args()

strings2check = ["DP:AD","AD:DP"]

list_of_folders = glob.glob(args.input+'*/')

for folder in list_of_folders:
	StrainName = folder.split('/')[-2]
	input_file = folder + StrainName + '_SNP/'+StrainName+'_variants.vcf'
	if os.path.isfile(input_file):
		# read vcf file from bcftools isec
		with open(input_file,"r")as f: 
			reader = csv.reader(f, delimiter="\t")
			vcf_in = list(reader)
		f.close()
		
		# remove headers etc
		vcf = []
		
		for i in vcf_in:
			#print i[0][0:2]
			if i[0][0:2] != "##":
				vcf.append(i)
		
		new = []
		for snp in vcf[1:]:
			chro = snp[0]
			posi = snp[1]
			refe = snp[3]
			alte = snp[4]
			qual = snp[5]
			formatdata = snp[8]
			for i in strings2check:
				if i in formatdata:
					#print snp[9]
					if args.format == "4.2":
						#print snp[9].split(":")[1]
						total = float(snp[9].split(":")[1]) #Total read depth DP
						freq = snp[9].split(":")[2].split(",") #Frequency of Reference and alternative alleles
					else:
						total = float(snp[9].split(":")[2]) #Total read depth DP
						if total == 0.0:
							total = 0.000000000000001
						freq = snp[9].split(":")[1].split(",")
						#print freq
						#print total
					if len(freq) == 2 :
						freqR = float(freq[0])/total
						freqA = float(freq[1])/total
						linea = [chro,posi,refe,alte,qual,total,freqR,freqA]
						new.append(linea)
					elif len(freq) == 3:
						# divide en dos
						alte1 = alte.split(",")[0]
						freqR1 = float(freq[0])/total
						freqA1 = float(freq[1])/total
						linea1 = [chro,posi,refe,alte1,qual,total,freqR1,freqA1]
			
						alte2 = alte.split(",")[1]
						freqR2 = float(freq[0])/total
						freqA2 = float(freq[2])/total
						linea2 = [chro,posi,refe,alte2,qual,total,freqR2,freqA2]
			
						new.append(linea1)
						new.append(linea2)
			
					else :
						print snp	
				
		with open(StrainName+".qaf.tab", "wb")as ofile :
			writer = csv.writer(ofile, delimiter = "\t")
			writer.writerow(["CHR","POS", "REF","ALT","QUAL","TOTC","AFR","AFA"])
			for row in new :
				writer.writerow(row)		
		
		#Change the PATH of the script
		snp_dist_cm1 = "Rscript ~/software/scripts/snp_distrib_v2.r " + StrainName+".qaf.tab " + args.chromos
		os.system(snp_dist_cm1)
		print StrainName + ' was done!'

print "DONE!"