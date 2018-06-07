#!/usr/bin/python

fname = "H37Rv_genome.fasta"
fh = open(fname, 'r')
fh.readline() # remove the fasta header
numA = 0
numT = 0
numC = 0
numG = 0
for line in fh:
    numA += line.count("A")
    numT += line.count("T")
    numC += line.count("C")
    numG += line.count("G")
allbases = float(numA+numT+numC+numG)
print ("sum:", numA+numT+numC+numG)
print ("A:", numA, "\tT:", numT, "\tC:", numC, "\tG:", numG)
print ("A:", numA/allbases, "\tT:", numT/allbases, "\tC:", numC/allbases, "\tG:", numG/allbases)

