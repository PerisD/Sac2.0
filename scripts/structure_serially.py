#! /usr/bin/env python
#coding: utf-8

##########################################################################################
#Script for STRUCTURE serially.                                                          #
#Authors: David Peris UW-Madison, Dept Genetics                                          #
#Usage: python structure_serially.py KGROUPS REPEATS OUTPUTFOLDER MAINPARAMS EXTRAPARAMS #
#                                                                                        #
#KGROUPS the number of groups to perform (MIN-MAX, e.g 1-8)                              #
#REPEATS the number of runs per groups                                                   #
#MAINPARAMS and EXTRAPARAMS are the files with the run parameters                        #
##########################################################################################

import sys,os
import argparse

helptext="""
This script is to perform Structure analysis serially
Authors: David Peris UW-Madison, Dept Genetics
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-k","--kgroups", help="number of K groups to be tested, example 1-8", type = str, default = '1-8')
parser.add_argument("-r","--repeats", help="number of repeats for each K", type = int, default = 5)
parser.add_argument("-o","--output", help="the output folder", type = str, default = None)
parser.add_argument("-m","--mainparams", help="the mainparams file", type = str, default = "mainparams")
parser.add_argument("-e","--extraparams", help="the extraparams file", type = str, default = "extraparams")

parser.set_defaults(spades=True)

args = parser.parse_args()

if not os.path.exists(args.output):
	os.makedirs(args.output)

minK = int(args.kgroups.split('-')[0])
maxK = int(args.kgroups.split('-')[1])

while minK <= maxK:
	counter = 1
	K = str(minK)
	minK += 1
	print K
	while counter <= args.repeats:
		print counter
		structure_cm1 = 'structure -K ' + str(K) + ' -o ' + args.output + 'K' + str(K) + '_' + str(counter)
		structure_cm1 += ' -m ' + args.mainparams + ' -e ' + args.extraparams + ' > ' + args.output + 'K' + str(K) + '_' + str(counter) + '.log'
		print 'Structure cm:'+structure_cm1
		os.system(structure_cm1)
		counter += 1

#ADD A STEP TO ZIP THE OUTPUTFILES FOR STRUCTURE HARVESTER

print "done" 
