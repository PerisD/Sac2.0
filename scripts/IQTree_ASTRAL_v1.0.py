__author__ = 'Peris'

#! /usr/bin/env python
#coding: utf-8

helptext="""
This script is to serially run IQTree for ASTRAL analysis
Authors: David Peris University of Oslo, Dept Biosciences
"""

import os
import argparse

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="an input text file with fasta genes to parse", type = str, default = None)
parser.add_argument("-f","--folder", help="folder where fasta are stored", type = str, default = None)
parser.add_argument("-e","--EndName", help="the last part of the name of the fasta file, i.e. _nt_tr.fas", type = str, default = None)
parser.add_argument("-t","--Threads", help="number of threads to be used, default 8", type = str, default = "2")
parser.add_argument("-s","--seed", help="seed number to run iqtree", type = str, default = "225494")

args = parser.parse_args()

if not os.path.exists("SLURM-outputs/"):
	os.makedirs("SLURM-outputs/")

fasta_list = open(args.input, "r")
fasta_list.next() #Skip the first line because it contains a header

for genes in fasta_list:
	gene_name = genes.split(',')[0]
	gene_name += args.EndName
	input_fasta = args.folder + gene_name
	iqtree_cmd = "iqtree -s " + input_fasta + " -B 1000 --wbt -T " + args.Threads + " --seed " + args.seed
	print iqtree_cmd
	os.system(iqtree_cmd)

print "Done!"
