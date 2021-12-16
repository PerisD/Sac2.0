#! /usr/bin/env python
#coding: utf-8

import sys,os
import argparse

helptext="""
This script will run AAF
Authors: David Peris UW-Madison, Dept Genetics
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="PATH to the folder where library folders are stored", type = str, default = None)
parser.add_argument("-o","--output", help="Outputname", type = str, default = None)
parser.add_argument("-k","--kmer", help="Number of kmers to be used", type = str, default ="21")

args = parser.parse_args()
parser.set_defaults(output=True)

input_folder = args.input
run_name = args.output
kmer = args.kmer

def aaf_phylokmer(input_folder,kmer,run_name): 
	cmdA = 'python2.7 ~/software/AAF/aaf_phylokmer.py -k ' + kmer + ' -t 4 -G 64 -n 3 '
	cmdA += '-d ' + input_folder + ' -o ' + run_name
	return cmdA

def aaf_distance(run_name): 
	cmdB = 'python2.7 ~/software/AAF/aaf_distance.py -i ' + run_name + '.gz'
	cmdB += ' -t 8 -G 1 -f ' + run_name + '.wc '
	cmdB += '-o ' + run_name
	return cmdB

os.system(aaf_phylokmer(input_folder,kmer,run_name))
os.system(aaf_distance(run_name))

print "Done!"