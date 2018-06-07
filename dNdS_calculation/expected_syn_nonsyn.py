#!/usr/bin/python

from genetic_code import code

# start of extractSeq ##########################
def extractSeq( fname ): # extract gene or protein seq from respective file handles
    fh = open(fname, 'r')
    seqdict = {}
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
                gname = gname+"_c" # not all complement genes end with a 'c'. Thus add this distinguishing feature to the gene name
            elif locstr.startswith("order"): # only one gene has this. [locus_tag=Rv3216] [location=order(3593369..3593437,3593439..3593852)]
                # enter this gene manually, encompassing the entire region covered by this gene
                locstr = "3593369..3593852"
            start, end = locstr.split('..')
            start = start.lstrip('<') # some gene start positions have location mentioned as [location=<2550340..2551326]. Remove the "<" sign.
            end = end.lstrip('>') # some gene end positions have location mentioned as 817531..>817866. Remove the '>' sign.
            #print(gname, start, end)
            seq = ""
        else:
            seq += line
    return seqdict
# end of extractSeq ##########################

# start of getH37Rv_proteins #################
def getH37Rv_proteins(fname):
    fh = open(fname, 'r')
    protnames = [] # names of protein coding genes only (no tRNA genes)
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
                    protnames.append(gname+"_"+locustag+"_c")
                else:
                    protnames.append(gname+"_"+locustag)
                flag = 0
                gname = ''
                locustag = ''
    return protnames
# end of getH37Rv_proteins #################

# start of computeNS #########################
def computeNS():
    NS = {} # key: codon. value: list of nonsyn and syn sites in the codon
    nt = ["A", "T", "C", "G"]
    for c in code:
        aa = code[c]
        cbases = list(c) # bases in the codon
        #print(aa, cbases)
        N = 0.0 # num nonsyn sites in this codon
        S = 0.0 # num syn sites in this codon
        for i in range(len(cbases)):
            n_n = 0 # num of nonsyn changes at this pos in codon
            n_s = 0 # num of syn changes
            for n in nt:
                if cbases[i] == n:
                    continue
                else: # mutate this base in the codon with the three other bases
                    mutcodon = cbases[:]
                    mutcodon[i] = n 
                    #print("\t",mutcodon)
                    mutaa = code["".join(mutcodon)]
                    if mutaa == aa:
                        n_s += 1   # synonymous change
                    else:
                        n_n += 1   # non synonymous change
            N += n_n/3.0
        S = 3 - N
        #print(N,S)
        NS[c] = [N, S]
    return NS
# end of computeNS #########################

# start of computeRefExpected ##############
def computeRefExpected( geneseq, codon_N_S, protnames):
    N_ref = 0.0 # expected nonsyn sites in gene sequences in H37Rv
    S_ref = 0.0 # expected syn sites in gene sequences in H37Rv
    for g in geneseq:
        if g in protnames: # only if this gene is a protein coding gene and not RNA or other type of gene
            seq = list(geneseq[g])
            #print(g)
            for i in range(3, len(seq), 3): # start from 2nd codon because many protein seq code for M but start with V, and then read gene sequence in steps of 3
                codon = "".join(seq[i:i+3])
                if codon in codon_N_S:
                    N_ref += codon_N_S[codon][0]
                    S_ref += codon_N_S[codon][1]
                else:
                    print("Error!! Codon doesn't exist.", codon, g)
    return([N_ref, S_ref])
# end of computeRefExpected ##############

geneseq = extractSeq( "H37Rv_genes.txt" ) # supply file name that has all gene seq in fasta format
protnames = getH37Rv_proteins( "H37Rv_genpept.gp" )
print("Num of genes:", len(geneseq))
print("Num of protein coding genes:", len(protnames))

codon_N_S = computeNS() # compute nonsyn and syn sites for each codon
genicregions_NS_H37Rv = computeRefExpected(geneseq, codon_N_S, protnames)
print (genicregions_NS_H37Rv)
