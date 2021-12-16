#!/usr/bin/env python
#coding: utf-8

import argparse
from subprocess import call
import os
import shutil

helptext="""
This script is a wrapper to go through the different steps performed by HybPiper pipeline
Authors: David Peris UW-Madison, Dept Genetics
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-t","--target", help="Fasta file with the target genes", type = str, default = None)
parser.add_argument("-i","--input", help="A tab tabulated text file with StrainName and path to the reads, remember to write * two read 1 and 2 PATH/StrainName_*.fastq o fq", type = str, default = None)
#parser.add_argument("-f","--format", help="Indicate if is nucleotide (dna) or aminoacid (aa)", type = str, default = "dna")
parser.add_argument("-y","--hybpiper", help="Indicate the path to the HybPiper folder", type = str, default = "~/software/HybPiper/")
parser.add_argument("-o","--output", help="Indicate the folder where the final multisequence fasta will be saved", type = str, default = ".")
parser.add_argument("-s","--sedY", type = str, default = "NO", help="There are cases where illumina reads includes in the line + the designation of pair 1 and 2 as +/1 or +/2 and it must to be removed, default NO")
parser.add_argument("-p","--spades", help="To avoid the Spades assembly, default YES",type = str, default = "YES")
parser.add_argument("-d","--depth", help="To perform coverage quantification, default NO",type = str, default = "NO")
parser.add_argument("--threads", type = str, default = "8", help="number of threads to be used, default 8")


parser.set_defaults()

args = parser.parse_args()
Strain_list = open("strain_list.txt","w")

if os.path.exists(args.output):
	print "The output directory %s already exists" % (args.output)
else:
	os.makedirs(args.output)

for line in open(args.input,"r"):
	StrainName = line.split("\t")[0]
	Reads = line.split("\t")[1].strip()
	if os.path.isfile(Reads.replace("*","1")):
		Strain_list.write(StrainName+"\n")
		if args.sedY == "YES":
			if not os.path.exists(args.output + '/modified_Iluminareads'):
				os.makedirs(args.output + '/modified_Iluminareads')
			old_Read1 = Reads.replace("*","1")
			new_Read1 = args.output + '/modified_Iluminareads/'+old_Read1.split('/')[-1]
			print new_Read1
			old_Read2 = Reads.replace("*","2")
			shutil.copy(old_Read1,new_Read1)
			sed_cmd = "sed -i 's/+\/1/+/g' " + new_Read1
			print sed_cmd
			os.system(sed_cmd)
			Reads = new_Read1.replace("_1","_*")
			Reads = os.getcwd()+"/"+Reads
			print Reads
			if os.path.isfile(old_Read2):
				new_Read2 = args.output + '/modified_Iluminareads/'+old_Read2.split('/')[-1]
				print new_Read2
				shutil.copy(old_Read2,new_Read2)
				sed_cmd = "sed -i 's/+\/2/+/g' " + new_Read2
				print sed_cmd
				os.system(sed_cmd)
		hybpiper_cmd = "python " + args.hybpiper + "reads_first.py -b " + args.target + " -r " + Reads + " --prefix " + StrainName 
		hybpiper_cmd += " --bwa --cpu " + args.threads + " --no-assemble --no-exonerate"
		print hybpiper_cmd
		os.system(hybpiper_cmd)
		if args.depth == "YES":
			if os.path.isfile(StrainName + "/" + StrainName + ".bam"):
				samtools_cmd1 = "samtools view --threads 7 -bS "+StrainName + "/" + StrainName+".bam > "+StrainName + "/" + StrainName+"_aln-pe.view.bam"
				print samtools_cmd1
				os.system(samtools_cmd1)
				samtools_cmd2 = "samtools sort "+StrainName + "/" + StrainName+"_aln-pe.view.bam -o "+StrainName + "/" + StrainName+"_aln-pe.sort.bam"
				print samtools_cmd2
				os.system(samtools_cmd2)
				geneCoverage_cmd = "/opt/bifxapps/bedtools2-2.27.0/genomeCoverageBed -d -ibam "+StrainName + "/" + StrainName+"_aln-pe.sort.bam > "+StrainName + "/" + StrainName+"_coverage.txt"
				print geneCoverage_cmd
				os.system(geneCoverage_cmd)

Strain_list.close()

get_seq_lengths_output = open("seq_lengths.txt","w")
hybpiper_stats_output = open("statistics_run.txt","w")

call(["python",args.hybpiper+"get_seq_lengths.py",args.target,"strain_list.txt","dna"],stdout=get_seq_lengths_output)
#In the future add the R script to show a heatmap
call(["python",args.hybpiper+"hybpiper_stats.py","seq_lengths.txt","strain_list.txt"],stdout=hybpiper_stats_output)
call(["python",args.hybpiper+"retrieve_sequences.py",args.target,args.output,"dna"])


if os.path.exists(args.output+"/read_hits"):
	print "The output directory %s already exists" % (args.output+"/read_hits")
else:
	os.makedirs(args.output+"/read_hits")

list_of_genes = []
genes_explored = open(args.target,"r")


for lines in genes_explored:
	if ">" in lines:
		gene_name = lines.split("-")[1].strip()
		list_of_genes.append(gene_name)

strain_explored = open("strain_list.txt","r")

for strain in strain_explored:
	strain = strain.strip()
	for gene in list_of_genes:
		if os.path.isfile(strain+"/"+gene+"/"+gene+"_interleaved.fasta"):
			shutil.copy(strain+"/"+gene+"/"+gene+"_interleaved.fasta",args.output+"/read_hits"+"/"+strain+"_"+gene+"_reads.fasta")
		elif os.path.isfile(strain+"/"+gene+"/"+gene+"_unpaired.fasta"):
			print "Moving %s's %s reads" % (strain,gene)
			shutil.copy(strain+"/"+gene+"/"+gene+"_unpaired.fasta",args.output+"/read_hits"+"/"+strain+"_"+gene+"_unpaired.fasta")
			print "Moved of %s's %s reads done!" % (strain,gene)

print "Done HybPiper Steps move to Spades!\n"

if args.spades == "YES":
	spades_directory = args.output+"/spades"
	if not os.path.exists(spades_directory):
		os.makedirs(spades_directory)
	for gene in list_of_genes:
		list4spades = open("strain_list.txt", "r")
		for strain in list4spades:
			strain_name = strain.strip()
			strain_directory = spades_directory + "/" + strain_name
			extracted_reads = args.output+"/read_hits/"+strain_name+"_"+gene+"_reads.fasta"
			call(["python","/opt/bifxapps/SPAdes-3.5.0/bin/spades.py","-t",args.threads,"-m","64","--12",extracted_reads,"-o",strain_directory,"--only-assembler"])
			if os.path.isfile(strain_directory+"/scaffolds.fasta"):
				shutil.copy(strain_directory+"/scaffolds.fasta",spades_directory+"/"+strain_name+"_"+gene+"_scaffolds.fas")

print "Spades step done!"
