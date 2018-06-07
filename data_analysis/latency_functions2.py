#!/usr/bin/python

import math
from genetic_code import code
import numpy as np, scipy.stats as st
#import matplotlib
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})

# define base complements
complement = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

##########################################
def readTruePairs( fname ): # read the file name that has the true pairs ids and the mutational distance between them
    truepairfh = open(fname, 'r')
    pairmutdist = {} # key: index case id, value: num mut diff between index and secondary case
    yeardict = {} # key: index case id, value: years between index and secondary case
    indexid = [] # list of index case ids
    secid = [] # corresponding secondary case ids
    flag = 0
    for line in truepairfh:
        if line.startswith("Idxcase"):
            flag = 1
            continue
        if flag == 1:
            line = line.rstrip("\n")
            content = line.split()
            try:
                idxid = int(content[0])
                mutdist = int(content[1])
                year = int(content[2])
                pairmutdist[idxid] = mutdist
                yeardict[idxid] = year
                indexid.append( idxid )
                secid.append( idxid+1 )
            except ValueError: # last line has result of the regression and is in float, not int
                break
    truepairfh.close()
    return pairmutdist, indexid, yeardict
##########################################

##########################################
def getMutData( infile, indexid ): # read csv file that has mutation pr/ab 1/0 for each SNP for each strain
    outfile = infile+".pair_mutations"
    fh = open(infile, 'r')
    fout = open(outfile, 'w')
    secid = [idx+1 for idx in indexid]

    header = fh.readline().split(',')
    # get column index matching the index-case and secondary-case ids
    idcolumn = {} # key: sample id (index or secondary), value: column id
    for i in range(1,len(header)-1): # 1st and last elts are "SNP" and "gname" so ignore those
        content = header[i].split('_')
        sid = int(content[0])
        if sid in indexid or sid in secid:
            idcolumn[sid] = i
        else:
            continue 

    #print (header)
    #for idx in indexid:
    #    print (idx, pairmutdist[idx], idcolumn[idx])
    #for idx in secid:
    #    print (idx, idcolumn[idx])

    # save mutations in each pair
    pairmutsnp = {} # key: index case id, value: list of mutations that are different in this index and its sec case
    pairmutgene = {} # key: index case id, value: list of genes that have mutations in this index and its sec case
    pairmutpos = {} # key: index case id, value: list of genome pos that have mutations in this pair 
    for idx in indexid: # initialize the above dicts
        pairmutsnp[idx] = []
        pairmutgene[idx] = []
        pairmutpos[idx] = []

    for line in fh:
        line = line.rstrip('\n')
        content = line.split(',')
        # get pairs where this snp differs in the index and sec cases
        for idx in indexid: # for each index case
            sid = idx+1 # secondary case id
            idxcol = idcolumn[idx]
            sidcol = idcolumn[sid]
            if content[idxcol] == content[sidcol]:
                continue # mutation either present or absent in both samples of a pair
            else:
                pairmutsnp[idx].append(content[0])
                pairmutgene[idx].append(content[-1])
                pairmutpos[idx].append(int(content[0][1:-1]))

    # print mutations in each pair
    for idx in indexid:
        fout.write("\nPair "+str(idx)+"-"+str(idx+1)+"\tMut dist:"+str(len(pairmutsnp[idx]))+"\n") #, set(pairmutgene[idx]))
        for i in range(len(pairmutsnp[idx])):
            fout.write(pairmutsnp[idx][i]+"\t"+pairmutgene[idx][i]+"\n")
    fout.close()
    return pairmutsnp, pairmutgene, pairmutpos
##########################################

# start of extractSeq ##########################
def extractSeq( fname ): # extract gene or protein seq from respective file handles
    fh = open(fname, 'r')
    seqdict = {}
    gposdict = {}
    seq = ""
    for line in fh:
        if ">" in line: # header line
            if seq != "": # there is some sequence to save
                seqdict[gname] = seq.replace('\n','')
            content = line.split()
            info = {'gene' : '', 'locus_tag' : '', 'location' : ''}
            for i in range(1, len(content)):
                key, val = content[i].split('=')
                key = key.lstrip('[')
                val = val.rstrip(']')
                if key in info: # H37Rv_proteins has a lot more fields, only get these three
                    info[key] = val
                else:
                    continue
            gname = info['gene']+"_"+info['locus_tag']
            locstr = info['location']
            if locstr.startswith("complement"): # info[2] is location string [location=12468..13016] or [location=complement(13133..13558)]
                locstr = locstr.lstrip("complement(")
                locstr = locstr.rstrip(")]")
                #gname = gname+"_c" # not all complement genes end with a 'c'. Thus add this distinguishing feature to the gene name
            elif locstr.startswith("order"): # only one gene has this. [locus_tag=Rv3216] [location=order(3593369..3593437,3593439..3593852)]
                # enter this gene manually, encompassing the entire region covered by this gene
                locstr = "3593369..3593852"
            start, end = locstr.split('..')
            start = start.lstrip('<') # some gene start positions have location mentioned as [location=<2550340..2551326]. Remove the "<" sign.
            end = end.lstrip('>') # some gene end positions have location mentioned as 817531..>817866. Remove the '>' sign.
            gposdict[gname] = [int(start), int(end)]
            #print(gname, start, end)
            seq = ""
        else:
            seq += line
    return seqdict, gposdict
# end of extractSeq ##########################

# start of getH37Rv_proteins #################
def getH37Rv_proteins(fname):
    fh = open(fname, 'r')
    protnames = [] # names of protein coding genes only (no tRNA genes)
    compgenes = [] # complement genes
    flag = 1
    gname = ''
    locustag = ''
    for line in fh:
        if line.startswith("LOCUS"):
            flag = 1
        if flag == 1:
            if "/gene=" in line:
                content = line.split('="')
                gname = content[1].rstrip('"\n')
            elif "/locus_tag" in line:
                content = line.split('="')
                locustag = content[1].rstrip('"\n')
            elif "/coded_by" in line:
                if "complement" in line: # complement gene
                    compgenes.append(gname+"_"+locustag)
                protnames.append(gname+"_"+locustag)
                flag = 0
                gname = ''
                locustag = ''
    return protnames, compgenes
# end of getH37Rv_proteins #################

# start of translateGene ##########################
def translateGene( nt ):
    tseq = ["M"] # no matter the first codon, the amino acid seq always starts with M -> NOT TRUE! Few proteins in H37Rv don't start with M
    for i in range(3, len(nt), 3): # start from 2nd codon, and then read gene sequence in steps of 3
        codon = "".join(nt[i:i+3])
        if codon in code:
            tseq.append( code[codon] )
        else: # some indels change frame, so last codon might be incomplete
            tseq.append( "***" )
    trans_aa = "".join(tseq)
    trans_aa = trans_aa.rstrip('*') # remove any stop codons at the end
    return trans_aa
# end of translateGene ##########################

##########################################
def getMutEffect( gname, gseq, snp, gposdict, compflag ):
    ntseq = list(gseq)
    pos = int(snp[1:-1])
    refbase = snp[:1]
    mutbase = snp[-1:]
    wtprot = translateGene(ntseq)
    gstart, gend = gposdict[gname]
    mutsite = pos - gstart
    mutseq = ntseq[:]
    if compflag == 1: # complement genes
        mutsite = -mutsite-1
        mutbase = complement[mutbase]
        refbase = complement[refbase]
        region = "".join(ntseq[mutsite-2:mutsite+2+1])
        #print (snp, gname, region)
    if ntseq[mutsite] == refbase:
        mutseq[mutsite] = mutbase
    else:
        print("ERROR! Refbase doesn't match!", snp, gname, "refbase_snp:", refbase, "refbase_seq:", ntseq[mutsite], "\n")
        print(ntseq[mutsite-2:mutsite+2+1])
    mutprot = translateGene(mutseq)
    wtseq = list(wtprot)
    mutatedseq = list(mutprot)
    aamut = ""
    for i in range(len(wtseq)):
        if wtseq[i] != mutatedseq[i]: 
            #print (i, wtseq[i], mutatedseq[i])
            aamut = wtseq[i]+str(i+1)+mutatedseq[i]
    if wtprot != mutprot:
        #print( gname, snp, "Non syn mutation:\n", wtprot, "\n", mutprot)
        return "nonsyn_"+aamut
    else:
        return "syn"
##########################################

# start of checkRegion ##########################
def checkRegion( pos, repeats ): # check if this position is in the repeat regions
    rregion = []
    for tup in repeats:
        if pos >= tup[0] and pos <= tup[1]: # pos is within the repeat region 
            rregion = tup
            break
    return rregion
# end of checkRegion ##########################

# start of compute_dNdS ########################
def compute_dNdS( Sd, Nd ): # compute dn/ds for each pair
    N = 2962658.6666573854 # reference expected nonsyn sites
    S = 1052142.333339442 # reference expected syn sites
    pN = Nd/N
    pS = Sd/S

    dN = -1 * 0.75 * math.log(1.0 - (4*pN)/3.0)
    dS = -1 * 0.75 * math.log(1.0 - (4*pS)/3.0)

    if dS == 0.0:
        return "Nan"
    else:
        return dN/dS 
# end of computedNdS ########################

# start of printMutArr #####################
def printMutArr(snplist):
    bases = ["A", "T", "C", "G"] # row and col headers in mutarr
    basemap = {"A" : 0, "T" : 1, "C" : 2, "G" : 3}
    mutarr = [ [0,0,0,0] for i in range(len(bases)) ]
    for snp in snplist:
        snparr = list(snp)
        wtbase = snparr[0]
        mutbase = snparr[-1]
        mutarr[basemap[wtbase]][basemap[mutbase]] += 1
    #print mutarr
    print("Row: WT, Col: Mutbase")
    print("  ", bases)
    nummut = 0
    for i in range(len(mutarr)):
        print(bases[i], ":", mutarr[i])
        nummut += sum(mutarr[i])
    print("Total num of mutations:", nummut)
# end of printMutArr #####################

# start of getOxidativeStressMut ####################
def getOxidativeStressMut(pairmutsnp, yeardict, indexid):
    year0to1_snps = []
    year2to5_snps = []
    print("========= Oxidative Stress Mutations =========")
    print(" Index ids: ", indexid)
    for idx in indexid:
        if idx == 49: # the pair that had 96 mutations
            continue
        year = yeardict[idx]
        if year == 0 or year == 1:
            year0to1_snps.extend(pairmutsnp[idx])
        else:
            year2to5_snps.extend(pairmutsnp[idx])
    print("Year 0 and 1 SNPs:")
    printMutArr(year0to1_snps)
    print("Year 2 to 5 SNPs:")
    printMutArr(year2to5_snps)
# end of getOxidativeStressMut ####################

# start of uniqueSNPsPerYear ################
def uniqueSNPsPerYear(pairmutsnp, yeardict, indexid, snpannotation):
    yearwisesnps = [[] for i in range(0,6)] # 2D list where index is year between pairs and value is all snps in pairs for that year.
    yearwisenewsnps = [[] for i in range(0,6)] # list where index is year between pairs and value is set of unique snps in pairs for that year.
    scatter_x = []
    scatter_y = []
    for idx in indexid:
        year = yeardict[idx]
        yearwisesnps[year].append(len(pairmutsnp[idx]))
        scatter_x.append(year)
        scatter_y.append(len(pairmutsnp[idx]))

    print ("SNPs each year:\n")
    avgsnp = [] # index: year difference between index and HHS, avg num SNPs for all pairs in that year
    boxplotdata = []
    for y in range(len(yearwisesnps)):
        avg_num_mut = float(sum(yearwisesnps[y]))/len(yearwisesnps[y])
        boxplotdata.append(yearwisesnps[y])
        avgsnp.append(avg_num_mut)
        print(y, yearwisesnps[y], avg_num_mut, "\n")
    #plt.plot(range(0,6),avgsnp, 'b')
    xticktext = ["$\leq$1\n$(n=7)$", "1-2\n$(n=4)$", "2-3\n$(n=3)$", "3-4\n$(n=6)$", "4-5\n$(n=4)$", "5-6\n$(n=1)$"]
    plt.boxplot(boxplotdata, 0, '+')
    plt.xticks(range(1,7),xticktext)
    #plt.axhline(y=0.0, color='gray', linestyle='--')
    plt.xlabel("Time between index case and HHC TB diagnosis (years)")
    plt.ylabel("Number of SNPs")
    plt.tight_layout()
    plt.savefig("Num_SNPs_boxplot.pdf") 
    plt.close()

    plt.plot(scatter_x, scatter_y, 'ko')
    plt.xticks(range(0,6),xticktext)
    #plt.axhline(y=0.0, color='gray', linestyle='--')
    plt.xlabel("Time between index case and HHC TB diagnosis (years)")
    plt.ylabel("Number of SNPs")
    plt.margins(0.05)
    plt.tight_layout()
    plt.savefig("Num_SNPs_scatterplot.pdf") 
    plt.close()
    return

# end of uniqueSNPsPerYear ################

# start of mean_ci_mu ###################
def mean_ci_mu(yeargroup, pairmutdist, timeinhours, gen, indexid, yeardict):
    N = 0.973 * 4411532 # 97.3% coverage for 4.41 MB genome of MTB. This is total number of sites in our analysis 
    mu_year01 = 1.75234300481e-08 #mu rate at 18 hr gen for year 0 and 1; old data: 1.91481189268e-08
    mu_range = mu_year01*2
    print (mu_year01-mu_range, mu_year01+mu_range, " <====")
    murates = {g : [] for g in gen} # mutation rates in the pairs that have year difference in list 'early'
    for idx in indexid:
        year = yeardict[idx]
        if year in yeargroup:
            print ("Year:", year, " index id:", idx)
            t = timeinhours[year]
            nsnps = pairmutdist[idx] 
            for g in gen:
                m = nsnps/(N*(t/g))
                murates[g].append(m) # append the mutation rate for this generation time

    yg_mean = []
    yg_95ci_upper = []
    yg_95ci_lower = []
    for g in gen:
        a = np.array(murates[g])
        a_mean = np.mean(a)
        if a_mean > mu_year01-mu_range and a_mean < mu_year01+mu_range: # identify generation that has active TB mu rate
            print (yeargroup, " mu_rate: ", a_mean, " gen: ", g)
        output = st.norm.interval(0.95, loc=a_mean, scale=np.std(a)/np.sqrt(len(a)))
        yg_mean.append(a_mean)
        yg_95ci_upper.append(output[1])
        yg_95ci_lower.append(output[0])
    return yg_mean, yg_95ci_upper, yg_95ci_lower
# end of mean_ci_mu ###################

# start of computeMuRate #################
def computeMuRate(pairmutdist, indexid, yeardict):
    N = 0.973 * 4411532 # 97.3% coverage for 4.41 MB genome of MTB. This is total number of sites in our analysis 
    gen = range(18,320,6) # generation time in hours. starts at 18 and ends at 240
    mu_year01 = 1.75234300481e-08 #mu rate for year 0 and 1; old data: 1.91481189268e-08
    early = [0,1] # early latency: years between index and sec case is 0 or 1
    late = [2, 3, 4, 5] # late latency
    year_grps = [[0,1],[2,3],[4,5]] # for mu-rate vs gen time figure
    timeinhours = {} # key: year, value: that much time in hours. Year 0 => 0.5 years => 0.5 x 365 x 24 = hours in 0.5 years
    for t in range(0,6):
        timeinhours[t] = float((t+0.5)) * 365 * 24

    f, axarr = plt.subplots(1, 3, sharey=True, figsize=(8,6))
    axarr[0].set_ylabel("Mutation rate\n(mutations/bp/generation)")
    for i in range(len(year_grps)):
        yg = year_grps[i]
        yg_mean, yg_95ci_upper, yg_95ci_lower = mean_ci_mu(yg, pairmutdist, timeinhours, gen, indexid, yeardict)
        print("year:", yg, " mu rate:", yg_mean[0])
        #print ("\n\nYear group: ", yg, "\nmean:", yg_mean, "\n95ci_upper:", yg_95ci_upper, "\n95ci_lower:", yg_95ci_lower)
        axarr[i].set_yscale('log')
        axarr[i].set_xlim([18,320])
        axarr[i].xaxis.set_ticks(np.arange(18, 320, 60))
        axarr[i].axhline(y=mu_year01, c='gray', linestyle='--')
        axarr[i].plot(gen, yg_mean, lw = 1, color = '#539caf', alpha = 1)
        axarr[i].fill_between(gen, yg_95ci_lower, yg_95ci_upper, color = '#539caf', alpha = 0.4)
        axarr[i].set_title("Years in latency: "+str(yg[0])+"-"+str(yg[1]+1), fontsize=10)
    f.text(0.5, 0.04, 'Generation time (hours)', ha='center')
    plt.savefig("latency_yeargroups_murates.pdf")
    #return

    gen18_murates = [] # mutation rates in each year group at gen=18 hours
    year_grps_yearwise = [[0], [1], [2], [3], [4], [5]] # for figures of mu rate decline and gen time increase, and linear regression
    gen_yearwise = range(18,750,6) # generation time in hours. starts at 18 and ends at 240
    for i in range(len(year_grps_yearwise)):
        yg2 = year_grps_yearwise[i]
        yg2_mean, yg2_95ci_upper, yg2_95ci_lower = mean_ci_mu(yg2, pairmutdist, timeinhours, gen_yearwise, indexid, yeardict)
        print (yg2, " gen 18 mean mu:", yg2_mean[0])
        gen18_murates.append(math.log(yg2_mean[0],10))
        print("year:", yg2, " mu rate:", yg2_mean[0])
    xvals = [0, 1, 2, 3, 4, 5] # yeargroups
    gen_at_wtmu = [18, 108, 126, 222, 408, 744] # generation time in each yeargroup that has mu rate closest to the active disease murate
    slope, intercept, r_val, p_val, serr = st.linregress(range(len(gen_at_wtmu)), gen_at_wtmu)
    plt.figure()
    plt.plot(range(len(gen_at_wtmu)),gen_at_wtmu, 'ko-')
    # plot the regression line
    x = range(len(gen_at_wtmu))
    y = [intercept + slope*val for val in x]
    #print("slope :", slope, " intercept: ", intercept)
    #print("+++ x:", x, "\n y :", y)
    plt.plot(x, y, c='gray', linestyle='--')
    plt.xlabel("Time between index case and HHC TB diagnosis (years)")
    plt.xticks(range(len(gen_at_wtmu)),xvals)
    plt.ylabel("Generation time (hours) corresponding to\nactive disease mutation rate")
    plt.margins(0.05)
    plt.tight_layout()
    plt.savefig("generation_time_vs_year.pdf")
    plt.close()
    print("Generation time vs latency years: linear regression slope and pval: ", slope, p_val)

    # plot the mutation rate vs latency years
    slope, intercept, r_val, p_val, serr = st.linregress(range(len(gen18_murates)), gen18_murates)
    plt.figure()
    plt.plot(range(len(gen18_murates)),gen18_murates, 'ko-')
    # plot the regression line
    x = range(len(gen18_murates))
    y = [intercept + slope*val for val in x]
    plt.plot(x, y, c='gray', linestyle='--')
    #plt.yscale('log')
    plt.xlabel("Time between index case and HHC TB diagnosis (years)")
    plt.xticks(range(len(gen18_murates)),xvals)
    plt.ylabel("Mutation rate ($log_{10}$ of mutations/bp/gen)\nat generation time = 18 hrs")
    plt.margins(0.05)
    plt.tight_layout()
    plt.savefig("mutation_rate_vs_year.pdf")
    plt.close()
    print("Mutation rate at gen=18hrs vs latency years: linear regression slope and pval: ", slope, p_val)

# end of computeMuRate ##################
