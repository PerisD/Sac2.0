__author__ = 'Quinn Langdon modified by Peris'

#!/usr/bin/env python
import sys

#Script to replace characters Y, R, K, M, S, W with a N in fasta

#fastaToChangePrefix = sys.argv[1]
#fastaToChangeName = fastaToChangePrefix + ".fasta"

fastaToChangeName = sys.argv[1]
outputFileName = sys.argv[2]

fastaToChange = open(fastaToChangeName, 'r')
lines = fastaToChange.readlines()
outFile = open(outputFileName, 'w')
for line in lines:
    currentLine = line.strip()
    if (">" in line) :
        outFile.write(currentLine + "\n")
    else:
        noY = currentLine.replace("Y", "N")
        noYR = noY.replace("R", "N")
        noYRK = noYR.replace("K", "N")
        noYRKM = noYRK.replace("M", "N")
        noYRKMS = noYRKM.replace("S", "N")
        noYRKMSV = noYRKMS.replace("V", "N")
        noAmbig = noYRKMSV.replace("W", "N")
        outFile.write(noAmbig + "\n")
outFile.close()