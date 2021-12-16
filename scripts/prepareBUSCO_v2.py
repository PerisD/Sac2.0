#! /usr/bin/env python
#coding: utf-8

import os
import shutil as sh
import argparse

helptext="""
This script will generate sbatch files for each assembly to run BUSCO v4 on them
Authors: David Peris UiO, Department of Biosciences
"""

parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", help="A text file with the STRAINNAME, ASSEMBLYPATH and species tag for augustus", type = str, default = None)
parser.add_argument("-t","--threads", help="Number of CPUs to be used", type = str, default = "8")
parser.add_argument("-l","--layoutFile", help="PATH to your BUSCO layout, choose the version of BUSCO", type = str, default = None)
parser.add_argument("-s","--dataBase", help="database to be used, make sure it was download first if not it \
will be downloaded several times, i.e. agaricomycetes_odb10", type = str, default = None)
parser.add_argument("-r","--repeatRun", help="If some runs failed and your are repeating them", type = str, default = "NO")
parser.add_argument("-d","--timeRun", help="time for the run, format: HH:MM:SS", type = str, default = "06:00:00")


parser.set_defaults()

args = parser.parse_args()

if not os.path.exists(os.getcwd()+"/SLURM-outputs/"):
	os.makedirs(os.getcwd()+"/SLURM-outputs/")

busco_bash = open("01_busco_bash.sh",'w')
list_assemblies2annotate = open(args.input)
noAssemblys = open("NoAssemblies_found.txt",'w')
for iAssembly in list_assemblies2annotate:
	StrainName = iAssembly.split('\t')[0]
	AssemblyPath = iAssembly.split('\t')[1]
	SppTAG = iAssembly.split('\t')[2].strip()
	if not os.path.isfile(AssemblyPath):
		noAssemblys.write(StrainName+"\t"+AssemblyPath+"\n")
	sh.copy(args.layoutFile,StrainName+"_busco.sh")
	busco_bash.write("sbatch " + StrainName+"_busco.sh\n")
	if args.repeatRun == "YES":
		os.system("rm -r " + StrainName + "/")
	modifySTRAINNAME = "sed -i 's|\[STRAINNAME\]|" + StrainName + "|g' " + StrainName+"_busco.sh"
	os.system(modifySTRAINNAME)
	modifyASSEMBLYPATH = "sed -i 's|\[ASSEMBLYPATH\]|" + AssemblyPath + "|g' " + StrainName+"_busco.sh"
	os.system(modifyASSEMBLYPATH)
	modifyTAG = "sed -i 's|\[SPPTAG\]|" + SppTAG + "|g' " + "*_busco.sh"
	os.system(modifyTAG)

modifyBUSCODB = "sed -i 's|\[BUSCODB\]|" + args.dataBase + "|g' " + "*_busco.sh"
os.system(modifyBUSCODB)
modifyTHREADS = "sed -i 's|\[THREADS\]|" + args.threads + "|g' " + "*_busco.sh"
os.system(modifyTHREADS)
modifyWORKINGDIR = "sed -i 's|\[WORKINGDIR\]|" + os.getcwd()+"/" + "|g' " + "*_busco.sh"
os.system(modifyWORKINGDIR)
modifyTIME = "sed -i 's|\[TIME\]|" + args.timeRun + "|g' " + "*_busco.sh"
os.system(modifyTIME)

busco_bash.close()
noAssemblys.close()
print("Done!")
