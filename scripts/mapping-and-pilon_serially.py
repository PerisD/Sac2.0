#!/usr/bin/env python
#coding: utf-8

import argparse
from subprocess import call
import os
import shutil

helptext="""
This script is to map the reads used in the assembly to correct the assembly using pilon, and to remap the reads to check for Heterozygosity values
Authors: David Peris UW-Madison, Dept Genetics
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A tab tabulated text file with StrainName and path to the genome assembly", type = str, default = None)
parser.add_argument("-t","--threads", help="Number of CPUs to be used", type = str, default = "8")
parser.add_argument("-q","--quality", help="Quality of mapping threshold", type = str, default = "30")
parser.add_argument("-e","--explore", help="If you only want to get quality of reads and assembly", type = str, default = "NO")


parser.set_defaults(spades=True)

args = parser.parse_args()

info_path = open(args.input, 'r')
info_path.readline()
wd_path = os.getcwd()

for line in info_path:
	line = line.split('\t')
	if line[2] == "YES":
		StrainName = line[0]
		print(StrainName)
		GenomePath = line[4]
		print(GenomePath)
		Read1 = line[5].strip()
		print(Read1)
		if '_1.fq' in Read1:
			Read2 = Read1.replace('_1.fq','_2.fq')
			print(Read2)
		if not os.path.exists(StrainName + '/'):
			os.makedirs(StrainName)
		if not os.path.isfile(StrainName + '/' + StrainName + '.fasta'):
			shutil.copy(GenomePath,StrainName + '/' + StrainName + '.fasta')
			print "Genome %s not found make a rename" % (StrainName)
		if not os.path.exists(StrainName + '/' + StrainName + '_SNP'):
			os.makedirs(StrainName + '/' + StrainName + '_SNP')
		print "=======================Starting pipeline %s========================" % (StrainName)
		bwa_cm1 = "bwa index -a is " + StrainName + '/' + StrainName + '.fasta'
		print bwa_cm1
		os.system(bwa_cm1)
		if '_1.fq' in Read1:
			bwa_cm2 = "bwa mem -t " + args.threads + " " + StrainName + '/' + StrainName + '.fasta ' + Read1 + " " + Read2 + " > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam"
		else:
			bwa_cm2 = "bwa mem -t " + args.threads + " " + StrainName + '/' + StrainName + '.fasta ' + Read1 + " > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam"
		print bwa_cm2
		os.system(bwa_cm2)
		samtools_cm1 = "samtools view -q " + args.quality + " -bhSu " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_view.sam"
		print samtools_cm1
		os.system(samtools_cm1)
		samtools_cm2 = "samtools sort -@ " + args.threads + " " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_view.sam -o " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam"
		print samtools_cm2
		os.system(samtools_cm2)
		samtools_cm3 = "samtools faidx ../" + StrainName + '/' + StrainName + '.fasta'
		print samtools_cm3
		os.system(samtools_cm3)
		if not os.path.exists(StrainName + '/' + StrainName + '_SNP/' + StrainName + '_qualimap'):
			os.makedirs(StrainName + '/' + StrainName + '_SNP/' + StrainName + '_qualimap')
		qualimap_cm1 = "~/software/qualimap_v2.2.1/qualimap bamqc -bam " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam -outdir " + StrainName + '/' + StrainName + '_SNP/' + StrainName + '_qualimap/'
		print qualimap_cm1
		os.system(qualimap_cm1)
		old_qual_name = StrainName + '/' + StrainName + '_SNP/' + StrainName + '_qualimap/genome_results.txt'
		print(old_qual_name)
		new_qual_name = StrainName + '/' + StrainName + '_SNP/' + StrainName + '_qualimap/'+ StrainName +'_genome_results.txt'
		print(new_qual_name)
		os.rename(old_qual_name,new_qual_name)
		if args.explore == "NO":
			picard_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/MarkDuplicates.jar I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam O="
			picard_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam M=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_picard-metrics.txt "
			picard_cm1 += "REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT"
			print picard_cm1
			os.system(picard_cm1)
			picard_cm2 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/AddOrReplaceReadGroups.jar  I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam O="
			picard_cm2 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_edup-ready.bam RGLB=runPEa RGPL=illumina RGSM=" + StrainName + " VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true RGPU=plateXXX"
			print picard_cm2
			os.system(picard_cm2)
			picard_cm3 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/CreateSequenceDictionary.jar R=" + StrainName + '/' + StrainName + '.fasta' +" O=" + StrainName + '/' + StrainName + ".dict"
			print picard_cm3
			os.system(picard_cm3)
			samtools_index = "samtools faidx " + StrainName + '/' + StrainName + ".fasta"
			print samtools_index
			os.system(samtools_index)
			gatk_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/gatk3/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + StrainName + '/' + StrainName + ".fasta -I "
			gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam --genotyping_mode DISCOVERY -mbq 20 -stand_emit_conf 31 -stand_call_conf 31 -o "
			gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf"
			print gatk_cm1
			os.system(gatk_cm1)
			pilon_cm1 = "/opt/bifxapps/jre7/bin/java -jar ~/software/pilon_v1.22/pilon-1.22.jar --diploid --iupac --changes --tracks --genome " + StrainName + '/' + StrainName + '.fasta ' + "--frags "
			pilon_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam --vcf --output " + StrainName + '/' + StrainName + '_SNP/' + StrainName + '_pilon'
			print pilon_cm1
			os.system(pilon_cm1)



print "DONE!"
