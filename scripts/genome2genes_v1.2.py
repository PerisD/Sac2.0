__author__ = 'Peris'

#!/usr/bin/env python
#coding: utf-8

import argparse
import os
import shutil
import glob
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import GC

helptext="""
This script is to extract genes and proteins from a genome assembly and aligned the orthologous genes based on table information
IMPORTANT!!!! Symbols like # is not allow in this version of the script
Authors: David Peris UW-Madison, Dept Genetics
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="Folder with genome assemblies and gff files", type = str, default = None)
parser.add_argument("-a","--annotation", help="A csv file indicating the mode to deal with the relationships ID SystName", type = str, default = None)
parser.add_argument("-c","--complete", help="When only one strain of the folder must be run but include the other strains in the final P/A gene list and alignments", type = str, default = "")


args = parser.parse_args()

list_of_genomes = glob.glob(args.input+'*.gff')

#Create lists of strains depending on complete ON or OFF
list_Strains2parse = []
list_of_genomes2parse = []
if args.complete == "":
	annotationFile = open(args.annotation,'r')
	annotationFile.next()
	for iLine in annotationFile:
		StrainName = Line.split(',')[0]
		list_of_genomes2parse.append(args.input+StrainName+'.gff')
		list_Strains2parse.append(StrainName)
	annotationFile.close()
elif args.complete != "":
	list_of_genomes2parse.append(args.input+args.complete+".gff")
	list_Strains2parse.append(args.complete)

print list_Strains2parse
print list_of_genomes2parse

###	Run process_gff_cds_proteins.pl from Jeffares

for gff_file in list_of_genomes2parse:
	StrainName = gff_file.split('/')[-1].split('.')[0]
	print("===========================Extracting genes from: "+StrainName+"===========================")
	NewFolder_indivgenes = StrainName + '_indivGenes'
	if not os.path.exists(NewFolder_indivgenes):
		os.makedirs(NewFolder_indivgenes)
	gff2genes_cmd1 = 'perl ~/software/scripts/process_gff_cds_proteins.pl -g ' + gff_file
	if os.path.isfile(args.input + StrainName+'.fasta'):
		gff2genes_cmd1 += ' -d ' + args.input + StrainName + '.fasta -o ' + NewFolder_indivgenes + '/'
	elif os.path.isfile(args.input + StrainName+'.fas'):
		gff2genes_cmd1 += ' -d ' + args.input + StrainName + '.fas -o ' + NewFolder_indivgenes + '/'
	print("gff2genes_cmd1:"+gff2genes_cmd1)
	os.system(gff2genes_cmd1)

###	Rename genes to be with identical systematic name

annotation_mode = pd.read_csv(args.annotation)
for Strain in list_Strains2parse:
	strain_info_temp = open(Strain+'_Info_renamed.txt','w')
	list_genes_rename = []
	list_genes_noRename = []
	Check_duplicates = []
	Indivgenes_folder = Strain + '_indivGenes/'
	Indivgenes_renamedFolder = Strain + '_indivGenes_renamed/'
	if not os.path.exists(Indivgenes_renamedFolder):
		os.makedirs(Indivgenes_renamedFolder)
	if len(os.listdir(Indivgenes_renamedFolder)) == 0:
		print Strain
		list_genes = glob.glob(Indivgenes_folder+'*dna.fasta') #FROM HERE REMOVE A TAB
		mode_anno = annotation_mode.loc[annotation_mode['StrainName'] == Strain, 'Annotation_Mode'].iloc[0]
		if mode_anno == "Peris":
			counterName = 1
			print "It is mode %s" % (mode_anno)
			annotation_gff = pd.read_csv(args.input+Strain+'.gff', sep ="\t", header = None)
			annotation_gff_temp = annotation_gff.loc[annotation_gff[2] == 'gene', 8].tolist()
			for gene in list_genes:
				GeneName = gene.split('/')[-1].split('.')[0]
				for item in annotation_gff_temp:
					if GeneName in item:
						if 'Script_Gene_selection' in item:
							item = item.split(';')
							for annotation in item:
								if 'Script_Gene_selection' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'Script_Gene_selection' in item and 'YGAP_Name' in item:
							item = item.split(';')
							for annotation in item:
								if 'YGAP_Name' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'YGAP_Name' in item and not 'Script_Gene_selection' in item:
							list_genes_noRename.append(GeneName)
		if mode_anno == "Peris2":
			counterName1 = 1
			counterName2 = 1
			counterName = 1
			print "It is mode %s" % (mode_anno)
			annotation_gff = pd.read_csv(args.input+Strain+'.gff', sep ="\t", header = None)
			annotation_gff_temp = annotation_gff.loc[annotation_gff[2] == 'gene', 8].tolist()
			for gene in list_genes:
				GeneName = gene.split('/')[-1].split('.')[0]
				for item in annotation_gff_temp:
					if GeneName in item:
						if '|' in item:
							item = item.split(';')
							for annotation in item:
								if '|' in annotation:
									SystName_temp = annotation.split('=')[1]
									GeneName1 = SystName_temp.split('|')[0]
									GeneName2 = SystName_temp.split('|')[1]
									if GeneName1 in Check_duplicates:
										Check_duplicates.append(GeneName1)
										GeneName1 = GeneName1+"_"+counterName1
										counterName1 += 1
									if GeneName2 in Check_duplicates:
										Check_duplicates.append(GeneName2)
										GeneName2 = GeneName2+"_"+counterName2
										counterName2 += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+GeneName1+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+GeneName1+'.prot.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+GeneName2+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+GeneName2+'.prot.fasta')
									list_genes_rename.append(GeneName1)
									list_genes_rename.append(GeneName2)
						elif not '|' in item and 'YGAP_Name' in item:
							item = item.split(';')
							for annotation in item:
								if 'YGAP_Name' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'YGAP_Name' in item and not 'Script_Gene_selection' in item:
							list_genes_noRename.append(GeneName)
		elif mode_anno == "Yue":
			counterName = 1
			print "It is mode %s" % (mode_anno)
			annotation_gff = pd.read_csv(args.input+Strain+'.gff', sep ="\t", header = None)
			annotation_gff_temp = annotation_gff.loc[annotation_gff[2] == 'gene', 8].tolist()
			for gene in list_genes:
				GeneName = gene.split('/')[-1].split('.')[0]
				for item in annotation_gff_temp:
					item = item.replace('_','',1)
					if GeneName in item:
						if 'Script_Gene_selection' in item:
							item = item.split(';')
							for annotation in item:
								if 'Script_Gene_selection' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'Script_Gene_selection' in item and 'Name' in item:
							item = item.split(';')
							for annotation in item:
								if 'Name' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'Name' in item and not 'Script_Gene_selection' in item:
							list_genes_noRename.append(GeneName)
		elif mode_anno == "Baker":
			counterName = 1
			print "It is mode %s" % (mode_anno)
			annotation_gff = pd.read_csv(args.input+Strain+'.gff', sep ="\t", header = None)
			#print annotation_gff #REMOVED ONCE DONE!
			annotation_gffGene_temp = annotation_gff.loc[annotation_gff[2] == 'gene', 8].tolist()
			annotation_gffCDS_temp = annotation_gff.loc[annotation_gff[2] == 'CDS', 8].tolist()
			#print annotation_gff_temp #REMOVED ONCE DONE!
			for gene in list_genes:
				GeneName = gene.split('/')[-1].split('.')[0]
				#print GeneName #REMOVED ONCE DONE!
				counterCDS = 0
				for itemGene in annotation_gffGene_temp:
					if GeneName in itemGene:
						#print GeneName #REMOVED ONCE DONE!
						if 'locus_tag' in itemGene:
							itemGene = itemGene.split(';')
							for annotationGene in itemGene:
								if 'locus_tag' in annotationGene:
									locus_tag = annotationGene.split('=')[1]
									#print locus_tag #REMOVED ONCE DONE!
									for itemCDS in annotation_gffCDS_temp:
										if locus_tag in itemCDS and counterCDS == 0:
											itemCDS = itemCDS.split(';')
											for annotationCDS in itemCDS:
												if 'saccharomyces cerevisiae ortholog' in annotationCDS:
													SystName_temp = annotationCDS.split('saccharomyces cerevisiae ortholog:')[1]
													SystName_temp = SystName_temp.replace(' ','')
													SystName_temp = SystName_temp.replace('"','')
													if 'hypotheticalprotein' in SystName_temp:
														SystName_temp = SystName_temp.replace('hypotheticalprotein','')
													#print SystName_temp #REMOVED ONCE DONE!
													counterCDS += 1
													if SystName_temp in Check_duplicates:
														Check_duplicates.append(SystName_temp)
														SystName_temp = SystName_temp+"_"+counterName
														counterName += 1
													shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
													shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
													list_genes_rename.append(GeneName)
										elif not 'saccharomyces cerevisiae ortholog' in itemCDS:
											list_genes_noRename.append(GeneName)
		elif mode_anno == "NCBI":
			counterName = 1
			print "It is mode %s" % (mode_anno)
			annotation_gff = pd.read_csv(args.input+Strain+'.gff', sep ="\t", header = None)
			annotation_gff_temp = annotation_gff.loc[annotation_gff[2] == 'gene', 8].tolist()
			for gene in list_genes:
				GeneName = gene.split('/')[-1].split('.')[0]
				for item in annotation_gff_temp:
					if GeneName in item:
						if 'Script_Gene_selection=' in item:
							item = item.split(';')
							for annotation in item:
								if 'Script_Gene_selection=' in annotation:
									SystName_temp = annotation.split('=')[1]
									if SystName_temp in Check_duplicates:
										Check_duplicates.append(SystName_temp)
										SystName_temp = SystName_temp+"_"+counterName
										counterName += 1
									shutil.copy(Indivgenes_folder+GeneName+'.dna.fasta',Indivgenes_renamedFolder+SystName_temp+'.dna.fasta')
									shutil.copy(Indivgenes_folder+GeneName+'.prot.fasta',Indivgenes_renamedFolder+SystName_temp+'.prot.fasta')
									list_genes_rename.append(GeneName)
						elif not 'Script_Gene_selection=' in item:
							list_genes_noRename.append(GeneName)
		strain_info_temp.write('The number of genes renamed for %s is %i\n'%(Strain,len(list_genes_rename)))
		for genes in list_genes_rename:
			genes = genes + '\n'
			strain_info_temp.write(genes)
		strain_info_temp.write('The number of genes no renamed is %i\n'%(len(list_genes_noRename)))
		for genes in list_genes_noRename:
			genes = genes + '\n'
			strain_info_temp.write(genes)
		strain_info_temp.write('The number of genes renamed more than 1 is %i\n'%(len(Check_duplicates)))
		for genes in Check_duplicates:
			genes = genes + '\n'
			strain_info_temp.write(genes)
		strain_info_temp.close()

### Renaming the sequence of each gene and protein

for StrainName in list_Strains2parse:
	print StrainName
	temp_path = StrainName + "_indivGenes_renamed/"
	list_genes = glob.glob(temp_path+'*dna.fasta')
	for file in list_genes:
		GeneName = file.split('/')[-1].split('.')[0]
		protFile = file.replace('dna','prot')
		fasta_file = open(file,'r')
		for line in fasta_file:
			if line.startswith('>'):
				old_line = line.strip()
				new_line = ">"+StrainName
		fasta_file.close()
		modify_dna = "sed -i 's/"+old_line+"/"+new_line+"/g' " + file
		modify_prot = "sed -i 's/"+old_line+"/"+new_line+"/g' " + protFile
		os.system(modify_dna)
		os.system(modify_prot)

list_Strains = []
annotationFile = open(args.annotation,'r')
annotationFile.next()
for iLine in annotationFile:
	StrainName = iLine.split(',')[0]
	list_Strains.append(StrainName)
annotationFile.close()

print list_Strains

###Making the total number of genes/proteins studied

list_totalGenes_Renamed = []
for StrainName in list_Strains:
	print StrainName
	temp_path = StrainName + "_indivGenes_renamed/"
	list_genes = glob.glob(temp_path+'*prot.fasta')
	for iFile in list_genes:
		GeneName = iFile.split('/')[-1].split('.')[0]
		if not GeneName in list_totalGenes_Renamed:
			if not "Umbiguous" in GeneName:
				list_totalGenes_Renamed.append(GeneName)

list_totalGenes_Renamed = sorted(list_totalGenes_Renamed)

table_PA = open('table_presence-absence.txt','w')
table_PA.write('GeneName\t'+'\t'.join(list_Strains)+'\n')

### Generating the file with P/A of genes

for GeneName in list_totalGenes_Renamed:
	list_presence = []
	list_presence.append(GeneName)
	for StrainName in list_Strains:
		if os.path.isfile(StrainName+"_indivGenes_renamed/"+GeneName+".prot.fasta"):
			value_temp = "1"
		else:
			value_temp = "0"
		list_presence.append(value_temp)
	table_PA.write('\t'.join(list_presence)+'\n')
table_PA.close()


print "DONE!"
