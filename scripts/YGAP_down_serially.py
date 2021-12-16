#!/usr/bin/env python
#coding: utf-8

import argparse
import os

helptext="""
This script is to download YGAP annotations projects
Authors: David Peris UW-Madison, Dept Genetics
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A tab tabulated text file with StrainName and name of the project", type = str, default = None)
parser.add_argument("-o","--output", help="output folder", type = str, default = None)

parser.set_defaults(spades=True)

args = parser.parse_args()

projectFile = open(args.input, 'r')

for line in projectFile:
	StrainName = line.split('\t')[0]
	ProjectName = line.split('\t')[1]
	projectNumber = line.split('\t')[2].strip()
	outputfolder = args.output + '/' + StrainName + '/'
	cmdYGAP1_Annotation = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/'+ProjectName+'.annotation_final.txt -P ' + outputfolder
	cmdYGAP2_tRNAannot = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/'+ProjectName+'.trna.txt -P ' + outputfolder
	cmdYGAP3_NoOrth = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/'+ProjectName+'.loneancestors.html -P ' + outputfolder
	cmdYGAP4_GetORF = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/getorf.html -P ' + outputfolder
	cmd_mvGetORF = 'mv ' + args.output + '/' + StrainName + '/getorf.html ' + args.output + '/' + StrainName + '/' + StrainName + '.getorf.html'
	cmdYGAP5_Dogs = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/dogs.html -P ' + outputfolder
	cmd_mvDogs = 'mv ' + args.output + '/' + StrainName + '/dogs.html ' + args.output + '/' + StrainName + '/' + StrainName + '.dogs.html'
	cmdYGAP6_Singletons = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/singletons.html -P ' + outputfolder
	cmd_Singletons = 'mv ' + args.output + '/' + StrainName + '/singletons.html ' + args.output + '/' + StrainName + '/' + StrainName + '.singletons.html'
	cmdYGAP7_GENBANK = 'wget http://wolfe.ucd.ie/annotation/finished/'+projectNumber+'/file/GENBANK.zip -P ' + outputfolder
	cmd_mvGENBANK = 'mv ' + args.output + '/' + StrainName + '/GENBANK.zip ' + args.output + '/' + StrainName + '/' + StrainName + '.GENBANK.zip'
	if not os.path.exists(args.output):
		os.makedirs(args.output)
	if not os.path.exists(outputfolder):
		os.makedirs(outputfolder)
	print "==============Download annotation data for %s" % (StrainName)
	os.system(cmdYGAP1_Annotation)
	os.system(cmdYGAP2_tRNAannot)
	os.system(cmdYGAP3_NoOrth)
	os.system(cmdYGAP4_GetORF)
	os.system(cmd_mvGetORF)
	os.system(cmdYGAP5_Dogs)
	os.system(cmd_mvDogs)
	os.system(cmdYGAP6_Singletons)
	os.system(cmd_Singletons)
	os.system(cmdYGAP7_GENBANK)
	os.system(cmd_mvGENBANK)

projectFile.close()

print "Done!"