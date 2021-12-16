__author__ = 'Peris'

#!/usr/bin/env python
#coding: utf-8

import argparse
from subprocess import call
import os
import shutil
import glob
from Bio import SeqIO

helptext="""
This script is to map reads to a reference assembly to generate an alignment for population genomics.
The output requires to convert N sites to - for posterior analysis using sed --i 's/N/-/g
Authors: David Peris UW-Madison, Dept Genetics & IATA-CSIC
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A tab tabulated text file with StrainName and read paths", type = str, default = None)
parser.add_argument("-t","--threads", help="Number of CPUs to be used, default 8", type = str, default = "8")
parser.add_argument("-q","--quality", help="Quality of mapping threshold, default 30", type = str, default = "30")
parser.add_argument("-r","--reference", help="Reference Genome path to map the reads", type = str, default = None)
parser.add_argument("-m","--maskText", help="Text File that will be used by maskCov_FASTA_byStrain-d_v2.1.py, default chromosome.txt", type = str, default = "chromosome.txt")
parser.add_argument("-W","--windowRHET", help="Window size to get the plots of Heterozygosity values, default is 10000 bp", type = str, default = "10000")
parser.add_argument("-p","--prefix", help="the prefix name added to the final alignment", type = str, default = "SNP")

parser.set_defaults(spades=True)

args = parser.parse_args()

GenomePath = args.reference

failed_strains = open("failed_strains.txt",'w')
failed_strains.write("Strains not parsed through the pipeline\n")

FOLDER2REFERENCE = GenomePath[:-len(GenomePath.split('/')[-1])]
Reference_StrainName = GenomePath.split('/')[-1].split('.')[0]

list_of_unambiguous_sequences = []

#Step to generate the text file that will be used in maskCov_FASTA_byStrain-d_v2.1.py
Chromosome_TextFile = open(args.maskText, "w")
for index, record in enumerate(SeqIO.parse(GenomePath,"fasta")):
	Chromosome_TextFile.write(record.id+"\n")
Chromosome_TextFile.close()

if not os.path.exists(GenomePath+'.amb'): #Generates the bwa index
	bwa_cm1 = "bwa index -a is " + GenomePath
	print "bwa_cm1:"+bwa_cm1
	os.system(bwa_cm1)
if not os.path.isfile(FOLDER2REFERENCE + Reference_StrainName + ".dict"): #Generates the picard dictionary
	picard_cm3 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/CreateSequenceDictionary.jar R=" + GenomePath +" O=" + FOLDER2REFERENCE + Reference_StrainName + ".dict"
	print "picard_cm3:"+picard_cm3
	os.system(picard_cm3)
if not os.path.isfile(GenomePath+'.fai'): #Generates the samtools index
	samtools_cm3 = "samtools faidx " + GenomePath
	print "samtools_cm3:"+samtools_cm3
	os.system(samtools_cm3)

info_path = open(args.input, 'r')
list_trains = []
for line in info_path:
	line = line.split('\t')[0]
	list_trains.append(line)
info_path.close()

list_trains_unique = []
for i in list_trains:
	if not i in list_trains_unique:
		globals()[i+"_counterLoop"] = 1
		globals()[i+"_accumulated_samtools_cm4"] = ""
		list_trains_unique.append(i)

for i in list_trains_unique:
	globals()[i] = list_trains.count(i)
	print i+" has:"+str(globals()[i])+" libraries"

info_path = open(args.input, 'r')
for line in info_path:
	line = line.split('\t')
	StrainName = line[0]
	Reads = line[1].split('\n')[0]
	Read1 = Reads.replace('*','1')
	Read2 = Reads.replace('*','2')
	#Section when just one library#
	if globals()[StrainName] == 1:
		print "=======================Starting pipeline %s========================" % (StrainName)
		print "=======================%s Library========================" % (globals()[StrainName])
		if not os.path.exists(StrainName + '/'):
			os.makedirs(StrainName)
		if not os.path.exists(StrainName + '/' + StrainName + '_SNP'):
			os.makedirs(StrainName + '/' + StrainName + '_SNP')
		if os.path.isfile(Read2):
			bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " " + Read2 + " > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam"
			print "bwa_cm2:"+bwa_cm2
			os.system(bwa_cm2)
		else:
			if os.path.isfile(Read1):
				bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam"
				print "bwa_cm2:"+bwa_cm2
				os.system(bwa_cm2)
			else:
				FOLDER2READ = Reads[:-len(Reads.split('/')[-1])]
				Read1 = glob.glob(FOLDER2READ+"*q")[0]
				bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam"
				print "bwa_cm2:"+bwa_cm2
				os.system(bwa_cm2)
		samtools_cm1 = "samtools view -q " + args.quality + " -bhSu " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".sam > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_view.sam"
		print "samtools_cm1:"+samtools_cm1
		os.system(samtools_cm1)
		samtools_cm2 = "samtools sort -@ " + args.threads + " " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_view.sam -o " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam"
		print "samtools_cm2:"+samtools_cm2
		os.system(samtools_cm2)
		picard_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/MarkDuplicates.jar I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam O="
		picard_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam M=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_picard-metrics.txt "
		picard_cm1 += "REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT"
		print "picard_cm1:"+picard_cm1
		os.system(picard_cm1)
		picard_cm2 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/AddOrReplaceReadGroups.jar  I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam O="
		picard_cm2 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam RGLB=runPEa RGPL=illumina RGSM=" + StrainName + " VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true RGPU=plateXXX"
		print "picard_cm2:"+picard_cm2
		os.system(picard_cm2)
		gatk_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/gatk3/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + GenomePath + " -I "
		gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam --genotyping_mode DISCOVERY -mbq 20 -stand_emit_conf 31 -stand_call_conf 31 -o "
		gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf"
		print "gatk_cm1:"+gatk_cm1
		os.system(gatk_cm1)
		if os.path.isfile(StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf"):
			VCF2FASTA_cm1 = "python ~/software/scripts/VCF-FASTAconvert.py " + GenomePath + " " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf " + StrainName + '/' + StrainName + "_" + args.prefix
			print "VCF2FASTA_cm1:"+VCF2FASTA_cm1
			os.system(VCF2FASTA_cm1)
			getHeterozygousSites_cm1 = "python ~/software/scripts/getHeterozygousSites-VCF.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf "
			getHeterozygousSites_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_HTZInfo"
			print "getHeterozygousSites_cm1:"+getHeterozygousSites_cm1
			os.system(getHeterozygousSites_cm1)
			plot_heterozygosity_cm1 = "Rscript ~/software/scripts/heterozygosityAverager+Plot.R " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_HTZInfo_genotype.txt "
			plot_heterozygosity_cm1 += args.windowRHET + " " + StrainName + '/' + StrainName + '_SNP/'
			print "plot_heterozygosity_cm1:"+plot_heterozygosity_cm1
			os.system(plot_heterozygosity_cm1)
			bedgraph_cm1 = "/opt/bifxapps/bedtools2-2.27.0/genomeCoverageBed -d -ibam " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".bedgraph"
			print "bedgraph_cm1:"+bedgraph_cm1
			os.system(bedgraph_cm1)
			DepthQuantiles_cm1 = "Rscript ~/software/scripts/depthQuantile_forMasking_byStrain_d.R " + StrainName + '/' + StrainName + '_SNP/' + StrainName
			print "DepthQuantiles_cm1:"+DepthQuantiles_cm1
			os.system(DepthQuantiles_cm1)
			MaskByCov_cm1 = "python ~/software/scripts/maskCov_FASTA_byStrain-d_v2.1.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".cov " + args.maskText + " "
			MaskByCov_cm1 += StrainName + '/' + StrainName + "_" + args.prefix + ".fasta " + GenomePath + ' pairwiseDivFile_chromosome.txt'
			print "MaskByCov_cm1:"+MaskByCov_cm1
			os.system(MaskByCov_cm1)
			unambiguous_cm1 = "python ~/software/scripts/FASTAremove_ambig_byStrain_v2.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_covMasked.fasta "
			unambiguous_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_unambiguous.fasta"
			print "unambiguous_cm1:"+unambiguous_cm1
			os.system(unambiguous_cm1)
			list_of_unambiguous_sequences.append(StrainName + '/' + StrainName + '_SNP/' + StrainName + "_unambiguous.fasta")
		else:
			failed_strains.write(StrainName+'\n')
	elif globals()[StrainName] > 1 and globals()[StrainName+"_counterLoop"] <= globals()[StrainName]:
		print "=======================%s Library========================" % (globals()[StrainName+"_counterLoop"])
		if not os.path.exists(StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/'):
			os.makedirs(StrainName + "_" + str(globals()[StrainName+"_counterLoop"]))
		if not os.path.exists(StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP'):
			os.makedirs(StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP')
		if os.path.isfile(Read2):
			bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " " + Read2 + " > " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + ".sam"
			print "bwa_cm2:"+bwa_cm2
			os.system(bwa_cm2)
		else:
			if os.path.isfile(Read1):
				bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " > " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + ".sam"
				print "bwa_cm2:"+bwa_cm2
				os.system(bwa_cm2)
			else:
				FOLDER2READ = Reads[:-len(Reads.split('/')[-1])]
				Read1 = glob.glob(FOLDER2READ+"*q")[0]
				bwa_cm2 = "bwa mem -t " + args.threads + " " + GenomePath + ' ' + Read1 + " > " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + ".sam"
				print "bwa_cm2:"+bwa_cm2
				os.system(bwa_cm2)
		samtools_cm1 = "samtools view -q " + args.quality + " -bhSu " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + ".sam > " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + "_view.sam"
		print "samtools_cm1:"+samtools_cm1
		os.system(samtools_cm1)
		samtools_cm2 = "samtools sort -@ " + args.threads + " " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + "_view.sam -o " + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + "_sort.bam"
		print "samtools_cm2:"+samtools_cm2
		os.system(samtools_cm2)
		globals()[StrainName+"_accumulated_samtools_cm4"] += StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + '_SNP/' + StrainName + "_" + str(globals()[StrainName+"_counterLoop"]) + "_sort.bam "
		print globals()[StrainName+"_accumulated_samtools_cm4"]
		globals()[StrainName+"_counterLoop"] += 1
	if globals()[StrainName] > 1 and globals()[StrainName+"_counterLoop"] > globals()[StrainName]:
		print "=======================Finally Merging Mappings %s Libraries========================" % (globals()[StrainName])
		print globals()[StrainName]
		print globals()[StrainName+"_counterLoop"]
		if not os.path.exists(StrainName + '/'):
			os.makedirs(StrainName)
		if not os.path.exists(StrainName + '/' + StrainName + '_SNP'):
			os.makedirs(StrainName + '/' + StrainName + '_SNP')
		samtools_cm4 = "samtools merge " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam " + str(globals()[StrainName+"_accumulated_samtools_cm4"])
		print "samtools_cm4:"+samtools_cm4
		os.system(samtools_cm4)
		picard_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/MarkDuplicates.jar I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_sort.bam O="
		picard_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam M=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_picard-metrics.txt "
		picard_cm1 += "REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT"
		print "picard_cm1:"+picard_cm1
		os.system(picard_cm1)
		picard_cm2 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/picard-tools-1.98/AddOrReplaceReadGroups.jar  I=" + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup.bam O="
		picard_cm2 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam RGLB=runPEa RGPL=illumina RGSM=" + StrainName + " VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true RGPU=plateXXX"
		print "picard_cm2:"+picard_cm2
		os.system(picard_cm2)
		gatk_cm1 = "/opt/bifxapps/jre7/bin/java -jar /opt/bifxapps/gatk3/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + GenomePath + " -I "
		gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam --genotyping_mode DISCOVERY -mbq 20 -stand_emit_conf 31 -stand_call_conf 31 -o "
		gatk_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf"
		print "gatk_cm1:"+gatk_cm1
		os.system(gatk_cm1)
		if os.path.isfile(StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf"):
			VCF2FASTA_cm1 = "python ~/software/scripts/VCF-FASTAconvert.py " + GenomePath + " " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf " + StrainName + '/' + StrainName + "_" + args.prefix
			print "VCF2FASTA_cm1:"+VCF2FASTA_cm1
			os.system(VCF2FASTA_cm1)
			getHeterozygousSites_cm1 = "python ~/software/scripts/getHeterozygousSites-VCF.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_variants.vcf "
			getHeterozygousSites_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_HTZInfo"
			print "getHeterozygousSites_cm1:"+getHeterozygousSites_cm1
			os.system(getHeterozygousSites_cm1)
			plot_heterozygosity_cm1 = "Rscript ~/software/scripts/heterozygosityAverager+Plot.R " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_HTZInfo_genotype.txt "
			plot_heterozygosity_cm1 += args.windowRHET + " " + StrainName + '/' + StrainName + '_SNP/'
			print "plot_heterozygosity_cm1:"+plot_heterozygosity_cm1
			os.system(plot_heterozygosity_cm1)
			bedgraph_cm1 = "/opt/bifxapps/bedtools2-2.27.0/genomeCoverageBed -d -ibam " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_dedup-ready.bam > " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".bedgraph"
			print "bedgraph_cm1:"+bedgraph_cm1
			os.system(bedgraph_cm1)
			DepthQuantiles_cm1 = "Rscript ~/software/scripts/depthQuantile_forMasking_byStrain_d.R " + StrainName + '/' + StrainName + '_SNP/' + StrainName
			print "DepthQuantiles_cm1:"+DepthQuantiles_cm1
			os.system(DepthQuantiles_cm1)
			MaskByCov_cm1 = "python ~/software/scripts/maskCov_FASTA_byStrain-d_v2.1.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + ".cov " + args.maskText + " "
			MaskByCov_cm1 += StrainName + '/' + StrainName + "_" + args.prefix + ".fasta " + GenomePath + ' pairwiseDivFile_chromosome.txt'
			print "MaskByCov_cm1:"+MaskByCov_cm1
			os.system(MaskByCov_cm1)
			unambiguous_cm1 = "python ~/software/scripts/FASTAremove_ambig_byStrain_v2.py " + StrainName + '/' + StrainName + '_SNP/' + StrainName + "_covMasked.fasta "
			unambiguous_cm1 += StrainName + '/' + StrainName + '_SNP/' + StrainName + "_unambiguous.fasta"
			print "unambiguous_cm1:"+unambiguous_cm1
			os.system(unambiguous_cm1)
			list_of_unambiguous_sequences.append(StrainName + '/' + StrainName + '_SNP/' + StrainName + "_unambiguous.fasta")
		else:
			failed_strains.write(StrainName+'\n')

failed_strains.close()

list_of_unambiguous_sequences.append(GenomePath)

outputTextFile_list = open("list_of_fastas.txt",'w')
for StrainFasta in list_of_unambiguous_sequences:
	outputTextFile_list.write(StrainFasta+"\n")
outputTextFile_list.close()

strain_chr_fastaUnambig_cm1 = "python ~/software/scripts/strain-chr-fastaUnambig.py list_of_fastas.txt " + args.maskText
print "strain_chr_fastaUnambig_cm1:"+strain_chr_fastaUnambig_cm1
os.system(strain_chr_fastaUnambig_cm1)

combineAllFASTA_cm1 = "python ~/software/scripts/combineAllFASTA_suffix.py list_of_fastas.txt AllGenomes"
print "combineAllFASTA_cm1:"+combineAllFASTA_cm1
os.system(combineAllFASTA_cm1)

print "DONE!"
