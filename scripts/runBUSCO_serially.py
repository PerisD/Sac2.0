#! /usr/bin/env python
#coding: utf-8

import argparse
import sys,os
from subprocess import call
import glob
import pandas as pd

helptext="""
This script will run BUSCO an generate a final report of all assemblies analyzed
Authors: David Peris UW-Madison, Dept Genetics
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-t","--threads", help="Number of CPUs to be used", type = str, default = "8")
parser.add_argument("-i","--input", help="PATH to the folder where assemblies are stored", type = str, default = None)
parser.add_argument("-o","--output", help="Folder to store the final runBUSCO output", type = str, default = "runBUSCO_report")

parser.set_defaults(output=True)

args = parser.parse_args()
assemblies_list = glob.glob(args.input+'*')

print assemblies_list

def run_BUSCO(assembly,StrainName):
	cmdA = 'python ~/software/busco/BUSCO.py -i ' + assembly 
	cmdA += ' -o ' + StrainName
	cmdA += ' -l ~/software/busco/saccharomycetales_odb9 -m geno'
	cmdA += ' -c ' + args.threads
	return cmdA

list_shortSummary = []
list_fullTable = []

for assembly in assemblies_list:
	print assembly
	StrainName = assembly.split('/')[-1].split('_')[0]
	if '.' in StrainName:
		StrainName = StrainName.split('.')[0]
	os.system(run_BUSCO(assembly,StrainName))
	list_shortSummary.append(os.getcwd()+'/run_'+StrainName+'/short_summary_'+StrainName+'.txt')
	list_fullTable.append(os.getcwd()+'/run_'+StrainName+'/full_table_'+StrainName+'.tsv')

counter = 1
for fullTable in list_fullTable:
	if counter == 1:
		StrainName = fullTable.split('/')[-1].split('_')[-1].split('.')[0]
		final_df_fullTable = pd.read_table(fullTable, skiprows =5, header = None, names = list('abcdefg'))
		final_df_fullTable.columns = ['Busco id',StrainName,'Contig','Start','End','Score','Length']
		final_df_fullTable = final_df_fullTable.ix[:,0:2]
		final_df_fullTable.to_pickle('./temp.txt')
		counter += 1
		print StrainName
	else:
		final_df_fullTable = pd.read_pickle('./temp.txt')
		StrainName = fullTable.split('/')[-1].split('_')[-1].split('.')[0]
		temporal_df = pd.read_table(fullTable, skiprows =5, header = None, names = list('abcdefg'))
		temporal_df.columns = ['Busco id',StrainName,'Contig','Start','End','Score','Length']
		final_df_fullTable = pd.merge(final_df_fullTable,temporal_df.ix[:,0:2], on = 'Busco id')
		final_df_fullTable.to_pickle('./temp.txt')
		print StrainName
		print counter

os.makedirs(args.output)
final_df_fullTable.to_csv(args.output+'/All_fullTable.txt', sep='\t')

final_df_shortSummary = pd.DataFrame(columns=('StrainName','Single-copy','Duplicated','Fragmented','Missing','Total searched'))
counter = 0
for shortSummary in list_shortSummary:
	StrainName = shortSummary.split('/')[-1].split('_')[-1].split('.')[0]
	temp_StrainName = StrainName
	temporal_shortSummary = open(shortSummary)
	for line in temporal_shortSummary:
		if 'single-copy' in line:
			temp_singlecopy = line.split('\t')[1]
		elif 'duplicated' in line:
			temp_duplicated = line.split('\t')[1]
		elif 'Fragmented' in line:
			temp_fragmented = line.split('\t')[1]
		elif 'Missing' in line:
			temp_missing = line.split('\t')[1]
		elif 'Total' in line:
			temp_total = line.split('\t')[1]
	temp_df = pd.DataFrame({'StrainName':[temp_StrainName],
                       'Single-copy':[temp_singlecopy],
                       'Duplicated':[temp_duplicated],
                       'Fragmented':[temp_fragmented],
                       'Missing':[temp_missing],
                       'Total searched':[temp_total]}, index = [counter])
	final_df_shortSummary = pd.concat([final_df_shortSummary,temp_df])
	counter += 1

final_df_shortSummary.to_csv(args.output+'/All_shortSummary.txt', sep='\t', index=False)

print "Done!" 
