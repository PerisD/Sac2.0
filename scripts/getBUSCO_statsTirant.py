#! /usr/bin/env python
#coding: utf-8

import argparse
import os
import pandas as pd
import shutil as sh

helptext="""
This script will get BUSCO data after running prepareBUSCO_v1 & sbatch
Authors: David Peris, Dept Biosciences, University of Oslo
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="list file with StrainName and PATH to assemblies used in prepareBUSCO_v1.py", type = str, default = None)
parser.add_argument("-f","--folderRun", help="full PATH to the folder where prepareBUSCO_v?.py was run", type = str, default = None)
parser.add_argument("-s","--dataBase", help="database used, i.e. agaricomycetes_odb10", type = str, default = None)
#parser.add_argument("-b","--buscoGenes", help="the file with information about genes in the used database, i.e. links_to_ODB10.txt", type = str, default = None)


parser.set_defaults()

args = parser.parse_args()
assemblies_list = open(args.input,'r')

list_shortSummary = []
list_fullTable = []

final_df_shortSummary = pd.DataFrame(columns=('StrainName','Single-copy','Duplicated','Fragmented','Missing','Total searched'))

if not os.path.exists("output_tables/"):
	os.makedirs("output_tables/")

counter1 = 1
for iStrain in assemblies_list:
	StrainName = iStrain.split('\t')[0]
	fullTable = args.folderRun+StrainName+'/run_'+args.dataBase+'/full_table.tsv'
	sh.copy(args.folderRun+StrainName+'/run_'+args.dataBase+'/full_table.tsv',"output_tables/"+StrainName+"_fullTable.tsv")
	shortSummary=args.folderRun+StrainName+'/run_'+args.dataBase+'/short_summary.txt'
	counter2 = 0
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
	temp_df = pd.DataFrame({'StrainName':[StrainName],
                       'Single-copy':[temp_singlecopy],
                       'Duplicated':[temp_duplicated],
                       'Fragmented':[temp_fragmented],
                       'Missing':[temp_missing],
                       'Total searched':[temp_total]}, index = [counter2])
	final_df_shortSummary = pd.concat([final_df_shortSummary,temp_df])
	counter2 += 1

final_df_shortSummary.to_csv('All_shortSummary.txt', sep='\t', index=False)

#combine_fullTable_cmd = "Rscript /storage/projects/vlc81/peris/software/scripts/Tirant/R_modules/busco_fullTable.R -i output_tables/ -b " + args.buscoGenes
#Using eukaryota_odb10 I do not see the links_to_ODB10.txt so I turn off this and the argument above and in the RScript
combine_fullTable_cmd = "Rscript /storage/projects/vlc81/peris/software/scripts/Tirant/R_modules/busco_fullTableTirant.R -i output_tables/"

print("Done!")
