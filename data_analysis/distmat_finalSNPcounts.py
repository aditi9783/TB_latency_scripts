#!/usr/bin/env python

# read the out.combined_highcov_SNPs_qscore200 files and compute SNP distance between all pairs of isolates. 
# 1st line is header that lists the order of strains.
# Each subsequent line starts with SNP and binary record of absence/presence of SNP in each of the strains. 
# Thus, find binary distance between each column. (Each col: binary vector of SNPs for a given strain) 

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 10})

fname = sys.argv[1] # file is "out.combined_highcov_SNPs_qscore200.allminus_PE_PPE"
fh = open(fname, 'r')
header = fh.readline().split(',')
#print header
sidx = [] # strains are not in order 1-60. Index of this list is the order they are currently in, and values in the list have order # from 1-60. (One pair is missing due to low data yield).
for i in range(1,len(header)-1): # 1st and last values are "SNP" and "gname". Ignore these.
    content = header[i].split('_')
    val = int(content[0])
    sidx.append(val)

#print sidx
binmat = []
for line in fh:
    content = line.split(',')
    binmat.append( [int(content[i]) for i in range(1,len(content)-1)] ) # 1st and last vals are SNP and genename, ignore those

binmatarr = np.array(binmat)
#print sidx[0], binmatarr[:,0] # print 1st column
#print sidx[1], binmatarr[:,1] # print 2nd column
#print sidx[-1], binmatarr[:,-1] # print last column
snpdict = {} # key: true sample id. Value: binary vector with snp presence/absence
for i in range(len(sidx)):
    snpdict[ sidx[i] ] = binmatarr[:,i]

# pair 3-4: RFLP dont match, 19-20: common strain, 45-46: non-consenting patient, 55-56: low coverage, 57-58: outlier with 776 SNPdiff
idrange = range(1,3)+range(5,19)+range(21,45)+range(47,55)+range(59,61) # sample ids are 1-60 but samples 55 and 56 are omitted because 56 had low data yield.
print idrange

###############################
# compute pairwise dist for all pairs
distmat = []
alldist = [] # list of all pairwise distances
oddnum = [idrange[i] for i in range(0,len(idrange)-1,2)] # the odd numbers in the sample list; the first of each paired samples
# get distances between paired samples as well as non pairs
paireddist = []
unpaireddist = []
for i in range(len(idrange)):
    distmat.append([0 for x in range(len(idrange))]) # initialize the distmat
for i in range(len(idrange)-1):
    s1 = idrange[i]
    for j in range(i+1,len(idrange)):
        s2 = idrange[j]
        d = np.count_nonzero(snpdict[s1]!=snpdict[s2])
        alldist.append(d)
        distmat[i][j] = d
        distmat[j][i] = d # symmetric matrix
        if s1 in oddnum and s2 == s1+1: # oddnum are 1st sample in pair, 2nd sample in pair is +1 of that. Eg, 1-2, 3-4, 5-6, .. are pairs
            paireddist.append(d)
        else: # not a pair
            unpaireddist.append(d)
            if d < 90:
                print s1, s2, d

#print "19v58", np.count_nonzero(snpdict[19]!=snpdict[58])
#print "19v35", np.count_nonzero(snpdict[19]!=snpdict[35])

#############################
# print distribution of distances
print "Number of true pairs:", len(paireddist), "Num all pairs minus the true pairs:", len(unpaireddist)
print "\nAverage SNP distance for true pairs: ", float(sum(paireddist))/len(paireddist)
print "Average SNP distance for non pairs: ", float(sum(unpaireddist))/len(unpaireddist)
hist_paired, bins_paired = np.histogram(paireddist, bins=50)
hist_unpaired, bins_unpaired = np.histogram(unpaireddist, bins=60)
print "True pairs:: histogram values:", hist_paired, "\n", "hist bins:", bins_paired, "\n\n"
print "Non true pairs:: histogram values:", hist_unpaired, "\n", "hist bins:", bins_unpaired
f = plt.figure()
plt.plot(bins_paired[0:-1], hist_paired, 'r', alpha=0.8)
plt.plot(bins_unpaired[0:-1], hist_unpaired, 'b', alpha=0.8)
plt.xlabel("SNP distance")
plt.ylabel("Frequency")
f.savefig("qscore200_SNPdist_No_outliers_pairedGreen_unpairedGray_hist.pdf")

###############################
# print distmat in phylip format
fout = open(fname+".with_57_58_pair_776SNPs_phylip_distmat", 'w')
numstrains = len(idrange)
fout.write("\t"+str(numstrains)+"\n")
counter = 1
for i in range(numstrains):
    sname = str(counter)+"i"
    if idrange[i] % 2 == 0: # divisible by zero, HHC case
        sname = str(counter)+"h"
        counter += 1
    # pad with spaces till strain name is 10 char long
    print sname
    while len(sname) < 10:
        sname += " "
    fout.write(sname+"\t"+"\t".join([str(v) for v in distmat[i]])+"\n")
fout.close()
