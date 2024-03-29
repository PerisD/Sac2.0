#!/usr/bin/env python

#Script generated by Jacek Kominek, UW-Madison, Dept Genetics

import os, sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='')
parser.add_argument("dir", help="Directory with BUSCO subdirs.")
parser.add_argument("--out", default="out", help="Output directory.")
parser.add_argument("--min", default="0", help="Directory with BUSCO subdirs.")
args = parser.parse_args()

busco = defaultdict(list)
total = 0
for l in sorted(os.listdir(args.dir)):
	if os.path.isdir(args.dir+"/"+l) == True and os.path.isdir(args.dir+"/"+l+"/single_copy_busco_sequences"):
		print l,
		ct = 0
		for f in os.listdir(args.dir+"/"+l+"/single_copy_busco_sequences"):
			if f.find("fna") != -1:
				seqs = list(SeqIO.parse(args.dir+"/"+l+"/single_copy_busco_sequences/"+f, "fasta"))
				if len(seqs) == 1:
					for s in seqs:
						#busco[f.split(".")[0]].append([f.split(".")[0]+l,str(s.seq)])
						busco[f.split(".")[0]].append([l,str(s.seq)])
						ct += 1
		print ct
		total += 1
print "Total",total

count_min = total
if int(args.min) > 0:
	count_min = int(args.min)

if os.path.exists(args.out) == False:
	os.makedirs(args.out)

for b in busco.keys():
	if len(busco[b]) >= count_min:
	#if len(busco[b]) == total:
		with open(args.out+"/"+str(len(busco[b]))+"_"+b+".fas","w") as busco_out:
			for f in busco[b]:
				busco_out.write(">"+f[0]+"\n"+f[1]+"\n")
	

		
