'''
Created on Jul 7, 2015
This is the new file for calculating the entropy and now onwards please use this. output is in the Entropy2015 on atlas
@author: rounakv
'''
from __future__ import division
import math, operator
import pysam, sys, pysamstats,csv 
import numpy as np
from Bio import SeqIO
import re


TaxonFile=csv.reader(open(sys.argv[4],"rU"), delimiter="\t")
# getting the genus from the taxonomy script
ContigGenus={} # stores contig as keys and genus as values
for row in TaxonFile:
    ContigGenus[row[0]]= re.split(" ",(re.split(";",row[1])[-1]))[0]

# calculate the average number of operations i.e mismatches, insertions or deletions normalized by coverage for that position 
def Operations(mismatCount,delCount,insertCount,depth):
    if mismatCount==0 and delCount==0 and insertCount==0:
        return(0)
    else:
        return((mismatCount+delCount+insertCount)/depth)
    
# Calculate Average Coverage and average operations over the contig length,  GC content,
samfile= pysam.Samfile(sys.argv[1],"rb")
fastafile=pysam.Fastafile(sys.argv[2]) 
ContigInfo={}
refs=samfile.references # contig ids
refLen=samfile.lengths # The lengths are in the same order as pysam.Samfile.references
rl=zip(refs,refLen)


# entropy file 
outfile=open(sys.argv[3],"wb")    
writefile=csv.writer(outfile,delimiter="\t")
writefile.writerow(["ContigName","ContigLength","AvgCoverage","GC%","AvgOperations","Genus"])

for ContigName, ContigLen in rl:
    print ContigName, ContigLen
    if ContigName not in ContigGenus:
        ContigGenus[ContigName]="unmapped"
    event=0
    coverage=0
    GC=0
    #print ContigName, ContigLen
    for StatsInfo in pysamstats.stat_variation(samfile,fastafile, chrom=ContigName, until_eof=True):
        event+=Operations(StatsInfo['mismatches'],StatsInfo['deletions'],StatsInfo['insertions'],StatsInfo['reads_all'])
        coverage+=StatsInfo['reads_all']
        if StatsInfo['ref']=="G" or StatsInfo['ref']=="C" or StatsInfo['ref']=="g" or StatsInfo['ref']=="c":
            GC+=1
            
    AvgOp=event/ContigLen # average no. of operations over the whole contig length
    AvgCov=coverage/ContigLen
    AvgGC=GC/ContigLen
    writefile.writerow([ContigName,ContigLen,AvgCov,AvgGC,AvgOp,ContigGenus[ContigName]])
outfile.close()

    
    