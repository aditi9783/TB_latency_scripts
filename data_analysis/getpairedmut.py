#!/usr/bin/python

from latency_functions2 import *

## MAIN ##
truepairf = "latency_pairs_mutdist_corrected_latency_dates.txt" # get file that has non-outlier pairs and num mut diff in each pair
pairmutdist, indexid, yeardict = readTruePairs(truepairf)

year_definition = {0:"<1", 1:"1-2", 2:"2-3", 3:"3-4", 4:"4-5", 5:"5-6"}

# read true mutation file and identify the true mutations different in each pair
infile = "out.combined_highcov_SNPs_qscore200.allminus_PE_PPE"
pairmutsnp, pairmutgene, pairmutpos = getMutData( infile, indexid )

allsnps = 0
for p in pairmutsnp:
    allsnps += len(pairmutsnp[p])
print ("All snps:", allsnps)

#ensembl_format_vep(pairmutsnp) # print SNP data in ensembl default format for running through the Variant Effect Predictor tool

geneseq, gposdict = extractSeq( "H37Rv_genes.txt" ) # supply file name that has all gene seq in fasta format
protnames, complementgenes = getH37Rv_proteins( "H37Rv_genpept.gp" ) # not all genes code for proteins. Get gene names of protein coding genes only.

#print (geneseq["_Rv0008c"])
#ntseq = list(geneseq["_Rv0008c"])
#trans = translateGene(list(geneseq["_Rv0008c"]))
#print (trans)
#protseq = readH37Rvproteins( "H37Rv_proteins.fasta" )

# test mutations in each pair
num_nonsyn = 0
num_syn = 0
nonsynsnps = [] # list of all non synonymous snps across all pairs
synsnps = []
snps_intergenic = [] # all intergenic snps
nonprotsnps = [] # snps in genes that do not code for proteins (for example tRNA genes)
pair_mutdist = {} # key: index id, val: list of syn, nonsyn, and intergenic snps
dNdS_dict = {} # key: index id, val: dN/dS
snpannotation = {} # key: each snp in all pairs, value: annotation as gene name or intergenic and syn or nonsyn mutation
nonsyn_mutgenes = {} # key: gene, value: list of index pairs that have nonsyn mut in that gene
syn_mutgenes = {} # key: gene, value: list of index pairs that have syn mut in that gene
for idx in indexid:
    print("\nPair "+str(idx)+"-"+str(idx+1)+"\tMut dist:"+str(len(pairmutsnp[idx]))+"\n") #, set(pairmutgene[idx]))
    this_nonsyn = 0
    this_syn = 0
    this_intergenic = 0
    this_nonprot_gene = 0 # genes that do not code for proteins (tRNA genes, for example)
    for i in range(len(pairmutsnp[idx])):
        snp = pairmutsnp[idx][i]
        gname = pairmutgene[idx][i]
        snpinfo = gname
        allgenes = []
        if ";" in gname: # overlapping genes. sometimes a mutation is in multiple genes because of that.
            allgenes = gname.split(';')
        else:
            allgenes = [gname]
        if gname != "intergenic":
            for g in allgenes:
                if g in protnames: # protein coding gene
                    gseq = geneseq[g]
                    cflag = 0 # cflag is 1 if this is a complement gene
                    if g in complementgenes:
                        cflag = 1
                    muttype = getMutEffect( g, gseq, snp, gposdict, cflag )
                    if muttype == "syn":
                        num_syn += 1
                        this_syn += 1
                        synsnps.append(snp)
                        if g in syn_mutgenes: # this gene has already found to have a nonsyn mutation
                            syn_mutgenes[g].append(idx)
                        else:
                            syn_mutgenes[g] = [idx]
                    elif muttype.startswith("nonsyn"):
                        num_nonsyn += 1
                        this_nonsyn += 1
                        nonsynsnps.append(snp)
                        if g in nonsyn_mutgenes: # this gene has already found to have a nonsyn mutation
                            nonsyn_mutgenes[g].append(idx)
                        else:
                            nonsyn_mutgenes[g] = [idx]
                    else:
                        print ("###### OUTLIER ##########")
                    print (snp, g, muttype)
                    snpinfo = snpinfo+"_"+muttype
                else: # gene is not protein coding
                    this_nonprot_gene += 1
                    nonprotsnps.append(snp)
        else:
            print (snp, gname, '#######')
            snps_intergenic.append(snp)
            this_intergenic += 1
        if snp not in snpannotation:
            snpannotation[snp] = snpinfo
        else:
            if snpannotation[snp] == snpinfo: # same snp found again
                continue
            else: # couple of snps were in diff genes, thus snpinfo should be different for each gene that has the same snp
                snpannotation[snp] = snpannotation[snp]+"_"+snpinfo
    pair_mutdist[idx] = [this_syn, this_nonsyn, this_nonprot_gene, this_intergenic]
    this_dNdS = compute_dNdS( this_syn, this_nonsyn )
    dNdS_dict[idx] = this_dNdS

print ("Pair\tYear\tLatencyYears\tTotal SNPs\tNum Synonymous\tNum Nonsynonymous\tNum RNA genes\tNum Intergenic\tdN/dS")
for idx in indexid:
    print(idx, "-", idx+1,"\t", yeardict[idx], "\t", year_definition[yeardict[idx]], "\t", sum(pair_mutdist[idx]), "\t", "\t".join([str(v) for v in pair_mutdist[idx]]), "\t", dNdS_dict[idx])

print ("Num syn:", len(synsnps), "\nNum nonsyn:", len(nonsynsnps), "\nNum non-protein coding genes:", len(nonprotsnps), "\nIntergenic:", len(snps_intergenic))

# print genes that have nonsyn mutations in multiple pairs
print ("\nGene Num_pairs_nonsyn_mut_in_gene")
singleton_nonsyn = [] # list of genes with a single nonsyn mutation in latency pairs
for g in nonsyn_mutgenes:
    numpairs_nonsynmut = len(set(nonsyn_mutgenes[g]))
    #if numpairs_nonsynmut > 1:
    print (g, numpairs_nonsynmut)
    #else:
    #    singleton_nonsyn.append(g)
print ("Total num genes with nonsyn mutations:", len(nonsyn_mutgenes))
print("Genes with single nonsyn mutation in pairs:", singleton_nonsyn)

# print genes that have syn mutations in multiple pairs
print ("\nGene Num_pairs_syn_mut_in_gene")
singleton_syn = [] # list of genes with a single syn mutation in latency pairs
for g in syn_mutgenes:
    numpairs_synmut = len(set(syn_mutgenes[g]))
    #if numpairs_synmut > 1:
    print (g, numpairs_synmut)
    #else:
    singleton_syn.append(g)
print ("Total num genes with syn mutations:", len(syn_mutgenes))
print("Genes with single syn mutation in pairs:", singleton_syn)

# print number of SNPs for each year minus the SNPs in year 0
uniqueSNPsPerYear(pairmutsnp, yeardict, indexid, snpannotation)

# make mutation rate plots as in Sarah Fortune's paper
computeMuRate(pairmutdist, indexid, yeardict)

# determine frequency of type of mutations in the early (year 0-1) and late (years 2-7) latency
getOxidativeStressMut(pairmutsnp, yeardict, indexid)

