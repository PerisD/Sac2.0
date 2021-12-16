__author__ = 'Peris'

#! /usr/bin/env python
#coding: utf-8

helptext="""
This script is to serially run IQTree
Authors: David Peris University of Oslo, Dept Biosciences
"""

import sys,os
import argparse
import glob
from Bio import SeqIO

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="an input folder where fas/fasta are stored, use the complete PATH or a text file with PATH to fasta to parse", type = str, default = None)
parser.add_argument("-b","--boots", help="Number of UltraFast bootstraps, choose 0 for non-bootstrap analysis, default 1000", type = str, default = "1000")
parser.add_argument("-t","--threads", help="Maximum number of threads to be used", type = str, default = "8")
parser.add_argument("-m","--model", help="Model to be used, default 'Test', IQTree will estimate the model", type = str, default = "Test")
parser.add_argument("-n","--network", help="You want the best trees in a Nexus file for Network/Cloudogram analysis (YES/NO), default NO", type = str, default = "NO")
parser.add_argument("-d","--trim", help="You want to trim your sequences using trimal NO == 0, proportion of sequences with gap,1 == no columns with gaps, or estimate\
by trimal with gappyout [0-1,gappyout] default 0", type = str, default = "0")
parser.add_argument("-s","--replacement", help="You need to replace gap symbol N by - (YES/NO), default NO", type = str, default = "NO")
parser.add_argument("-o","--outgroup", help="You want to root the tree using an outgroup", type = str, default = "")
#parser.add_argument("-f","--newfolder", help="Output folder", type = str, default = None)


args = parser.parse_args()

if args.input.endswith('/'):
	fasta_list = glob.glob(args.input+'*.fas*')
else:
	fasta_list = []
	temporal_list = open(args.input,'r')
	for iFasta in temporal_list:
		fasta_list.append(iFasta.strip())

#if not os.path.exists(args.newfolder):
#	os.makedirs(args.newfolder)
#output_dir = os.getcwd() + '/' + args.newfolder

#Trim sequences using trimal
if not args.trim == "0":
	new_fasta_list = []
	for genes in fasta_list:
		gene_name = genes.split('/')[-1]
		gene_name = gene_name.split('.')[0]
		genes = genes.strip()
		temp_folder = 'temp_folder'
		if args.replacement == "YES":
			if not os.path.exists(temp_folder):
				os.makedirs(temp_folder)
			fasta_seq = open(genes,'r')
			new_fasta_seq = open(temp_folder+'/'+gene_name+'.fas','w')
			for line in fasta_seq:
				if line.startswith('>'):
					new_fasta_seq.write(line)
				else:
					line = line.replace('N','-')
					new_fasta_seq.write(line)
			fasta_seq.close()
			new_fasta_seq.close()
			genes = temp_folder+'/'+gene_name+'.fas'
		if not args.trim == "gappyout":
			trimal_cmd1 = "trimal -in " + genes + " -out " + gene_name + "_tr.fas -gt " + args.trim
			print("trimal_cmd1:"+trimal_cmd1)
			os.system(trimal_cmd1)
		else:
			trimal_cmd1 = "trimal -in " + genes + " -out " + gene_name + "_tr.fas -gappyout"
			print("trimal_cmd1:"+trimal_cmd1)
			os.system(trimal_cmd1)
		new_fasta_list.append(gene_name + "_tr.fas")
	fasta_list = new_fasta_list

#IQTree analysis
if args.model == 'Test':
	model = 'MFP'
else:
	model = args.model

for genes in fasta_list:
	gene_name = genes.split('/')[-1].split('.')[0]
	gene_seq = genes.strip()
	iqtree_cmd1 = 'iqtree -nt AUTO -ntmax ' + args.threads + ' -B ' + args.boots + ' -alrt 1000 -m ' + model + ' -s ' + gene_seq + ' -pre ' + gene_name
	if not args.outgroup == "":
		iqtree_cmd1 += ' -o ' + args.outgroup
	print("iqtree_cmd1:"+iqtree_cmd1)
	os.system(iqtree_cmd1)

#Compile all best trees in one file
if args.network == "YES":
	output_file = "all_IQTrees.nex"
	newick_files = glob.glob('*.treefile')
	list_of_strains = []
	for genes in fasta_list:
		for strain_name in SeqIO.parse(genes, "fasta"):
			if not strain_name.id in list_of_strains:
				list_of_strains.append(strain_name.id)
	total_names = sorted(list_of_strains)
	print "Your total number of strains is %i" % (len(total_names))
	new_file = open(output_file, 'w')
	new_file.write("#nexus\n\n")
	new_file.write(("BEGIN Taxa;\nDIMENSIONS ntax=%i;\nTAXLABELS\n") % (len(total_names)))
	counter = 1
	for strains in total_names:
		new_file.write(("[%i] '%s'\n") % (counter, strains))
		counter += 1
	new_file.write(";\nEND; [Taxa]\n\n")
	new_file.write("BEGIN Trees;\nPROPERTIES partialtrees=yes;\nTRANSLATE\n")
	for strains in total_names:
		new_file.write(("\t '%s'\t'%s',\n") % (strains, strains))
	new_file.write(";\n")
	new_file.write("[TREES]\n")
	counter = 1
	for tree_file in newick_files:
		tree_name = tree_file.split('.')
		phylo_tree = (open(tree_file)).next()
		phylo_tree = phylo_tree[:-2]
		new_file.write(("[%i] tree '%s' = [&R] %s;\n") % (counter, tree_name[0], phylo_tree))
		counter = counter + 1
	print "The total number of trees is %i" % (counter-1)
	new_file.write("\nEND; [Trees]\n")
	new_file.close()
	print "New %s is ready to use" % (output_file)

print "Done!"
