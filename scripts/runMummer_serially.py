#! /usr/bin/env python
#coding: utf-8

#######################################################################################################
#Script for constructing dotplots using mummer.                                                       #
#Authors: David Peris UW-Madison, Dept Genetics                                                       #
#Usage: python runMummer.py ARG1                                                                      #
#                                                                                                     #
#ARG1: tabulated text file 1st column (Path2Assembly/Assembly name) (Path2Reference/Reference_Name)   #
#                                                                                                     #
#######################################################################################################

import sys,os
from subprocess import call

input_file = sys.argv[1]

list_of_assemblies = open(input_file)

def move_files(strain_name):
	cmdA = 'mv '+strain_name+'.* '+strain_name+'/'
	return cmdA

def show_coords(strain_name):
	cmdB = 'show-coords -r '+strain_name+'.delta'+' > '+strain_name+'.coords'
	return cmdB

for line in list_of_assemblies:
	line = line.strip()
	assembly_path = line.split('\t')[0]
	ref_path = line.split('\t')[1]
	strain_name = assembly_path.split('/')[-1].split('.')[0]
	new_folder = os.getcwd()+'/'+strain_name
	os.makedirs(new_folder)
	print " ".join(["nucmer","--maxgap=500","--mincluster=150","--prefix="+strain_name,ref_path,assembly_path])
	call(["nucmer","--maxgap=500","--mincluster=150","--prefix="+strain_name,ref_path,assembly_path])
	os.system(show_coords(strain_name))
	print " ".join(["mummerplot","--png","--prefix="+strain_name,"-R",ref_path,"-Q",assembly_path,"--filter",strain_name+'.delta'])
	call(["mummerplot","--png","--prefix="+strain_name,"-R",ref_path,"-Q",assembly_path,"--filter",strain_name+'.delta'])
	print " ".join(["mummerplot","--postscript","--prefix="+strain_name,"-R",ref_path,"-Q",assembly_path,"--filter",strain_name+'.delta'])
	call(["mummerplot","--postscript","--prefix="+strain_name,"-R",ref_path,"-Q",assembly_path,"--filter",strain_name+'.delta'])
	os.system(move_files(strain_name))

print "Done!" 
