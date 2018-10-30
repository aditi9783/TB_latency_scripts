#!/usr/bin/env python

# To determine which model best fits the observed data
# M1: No change in mutation accumulation. Rate param of poisson stays the same as year 0 (i.e. 14.2 SNPs)
# M2: Rate declines almost linearly with rate exp(-x/10) where x is years/
# M3: Rate decays exponentially with rate exp(-x)
# M4: Rate decays sharply exponentially with rate exp(-x*3)
# M5: Rate decays sharply exponentially with rate exp(-x*10)

import numpy as np, scipy.stats as st
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 12})


#avgsnp_yearwise = [14.75, 32.2, 16.5, 13.67, 13, 8.67, 17, 7] # this is from old dates
#N_yearwise = [4, 5, 4, 3, 4, 3, 1, 1] # from old dates

#avgsnp_yearwise = [26.28, 16.0, 14.75, 14.83, 11.33, 7.5] #new dates, including 96 SNPpair
#N_yearwise = [7, 3, 4, 6, 3, 2] #new dates, including 96 SNPpair

#avgsnp_monthwise = [14.5, 42, 17.5, 8.5, 24.5, 6, 15.33, 14.33, 15, 9, 7] # in 6-month scale
#N_monthwise = [4, 3, 2, 2, 2, 1, 3, 3, 1, 3, 1] # data in 6-month scale

# NEXT TWO LINES ARE FOR FINAL 24 pairs (MINUS VV3010)
avgsnp_monthwise = [14.5, 42, 17.5, 8.5, 24.5, 6, 15.33, 14.33, 15, 11.5, 7] # in 6-month scale
N_monthwise = [4, 3, 2, 2, 2, 1, 3, 3, 1, 2, 1] # data in 6-month scale

#print N_yearwise, N_yearwise_reactivation

N_yearwise_reverse = [4, 4, 4, 4, 5, 4] # num samples from which the avg values in avgsnp_yearwise_reverse was calculated
#avgsnp_yearwise = [14.67, 16.0, 14.75, 14.83, 11.33, 7.5] #new dates, omitting 96 SNPpair
#N_yearwise = [6, 3, 4, 6, 3, 2] #new dates, omitting 96 SNPpair

rate = avgsnp_monthwise[0] # this is the default rate of actively replicating bacteria i.e. rate in year 0

### start of getPoissonReactivation ##########
def getPoissonReactivation(rate): 
    x = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    y6 = []
    y7 = []
    plt.figure()
    for i in range(0, len(N_monthwise)):
        x.append(i)
        y1.append(rate * np.exp(0))
        y2.append(1.46 * float(i)) # linear. y=mx, for x=11 (last time point) y=14.5. Thus m ~ 1.32
        y4.append(1e-3 * np.exp(float(i) * 0.96))
        y5.append(1e-5 * np.exp(float(i) * 1.42))
        y6.append(1e-22 * np.exp(float(i) * 5.33))
    #plt.plot(x,y1,'r-', label="$\epsilon$ = 0")
    #plt.plot(x,y2,'g-', label="$\epsilon$ = -x/10")
    #plt.plot(x,y4,'b-', label="$\epsilon$ = -x")
    #plt.plot(x,y5,'k-', label="$\epsilon$ = -x*2")
    #plt.plot(x,y6,'y-', label="$\epsilon$ = -x*10")
    plt.plot(x,y1,'r-', label="Constant")
    plt.plot(x,y2,'g-', label="Linear")
    plt.plot(x,y4,'b-', label="Exponential")
    #plt.plot(x,y5,'k-', label="Rapid Exponential-1")
    plt.plot(x,y6,'y-', label="Rapid Exponential")
    #plt.legend(loc=1, bbox_to_anchor=(1.0, 0.95), fontsize=9)
    xticktext = ["5-5.5", "4.5-5", "4-4.5", "3.5-4", "3-3.5", "2.5-3", "2-2.5", "1.5-2", "1-1.5", "0.5-1", "$\leq$0.5"]
    #xticktext = ["$\leq$0.5","0.5-1", "1-1.5", "1.5-2", "2-2.5", "2.5-3", "3-3.5", "3.5-4", "4-4.5", "4.5-5", "5-5.5"]
    plt.xticks(range(0,len(N_monthwise)), xticktext, rotation=30) 
    plt.xlabel("Years before reactivation of TB in HHC")
    plt.ylabel("New SNPs in year")
    #plt.ylim([0,20])
    plt.tight_layout()
    #plt.title("Hypotheses of exponential decay of mutation incidence, $exp(\epsilon)$")
    plt.savefig("expgrowth_reactivation.pdf")
    plt.close()

    #print y1, "\n", y2, "\n", y4, "\n", y5, "\n", y6
    return [y1[::-1], y2[::-1], y4[::-1], y5[::-1], y6[::-1]] # because every reactivated patient has to go through reactivation phase, reverse the rates list to start drawing samples from reactivation and then walking back
### end of getPoissonGrowth ##########

### start of getPoissonGrowth ##########
def getPoissonGrowthBoth(avgsnps): 
    origbins= len(N_monthwise)
    #x = range(0,origbins*2)
    react_pt = range(1,22,2) # reactivation points for each bin of 6 months. Each bin furter modeled as two 3-month bins, thus 11x2 = 22 total time points with react_pt defining the reactivation points for each of the 6-month bins.
    binwise_rates_both = []
    for i in range(len(avgsnps)): # each value is avg snps in that 6-month bin.
        x = []
        rates = {"Constant" : [], "Linear" : [], "Expo" : [], "RapidExpo_2" : [], "RapidExpo_10" : []} # initialize rate arrays
        eachside_rate = float(avgsnps[i])/2.0 # assuming equal number of mutation are accrued during early latency and reactivation
        reactivation_pt = react_pt[i]
        totalbins = (i+1)*2 # avgsnps[0] has data for first 6months. totalbins = (0+1)x2 = 2 to get num total bins for this period in 3months scale
        #print i, "react pt", reactivation_pt, "es_rate:", eachside_rate
        for j in range(0, reactivation_pt+1):
        #    print j
            x.append(j)
            rates["Constant"].append(eachside_rate * np.exp(0))
            rates["Linear"].append(-1 * eachside_rate/(origbins*2) * float(j) + eachside_rate)
            rates["Expo"].append(eachside_rate * np.exp(-1 * j))
            rates["RapidExpo_2"].append(eachside_rate * np.exp(-1 * j * 2))
            rates["RapidExpo_10"].append(eachside_rate * np.exp(-1 * j * 10))
        #print i, rates["Linear"]
        #print "react pt:", rates["Linear"][reactivation_pt]
        #print "slope:", (eachside_rate-rates["Linear"][reactivation_pt])/(origbins*2 - reactivation_pt)
        for k in range(reactivation_pt+1, totalbins+1):
        #    print k
            x.append(k)
            rates["Constant"].append(eachside_rate * np.exp(0))
            intercept_Linear = rates["Linear"][reactivation_pt] - ((eachside_rate-rates["Linear"][reactivation_pt])/(totalbins - reactivation_pt))*reactivation_pt # y - mx, using reactivation_pt as the point from which y and x values are used
            rates["Linear"].append((eachside_rate-rates["Linear"][reactivation_pt])/(totalbins - reactivation_pt) * float(k) + intercept_Linear)
            intercept_Expo = rates["Expo"][reactivation_pt] - ((eachside_rate-rates["Expo"][reactivation_pt])/(totalbins - reactivation_pt))*reactivation_pt # y - mx, using reactivation_pt as the point from which y and x values are used
            rates["Expo"].append((eachside_rate-rates["Expo"][reactivation_pt])/(totalbins - reactivation_pt) * float(k) + intercept_Expo)
            intercept_RapidExpo_2 = rates["RapidExpo_2"][reactivation_pt] - ((eachside_rate-rates["RapidExpo_2"][reactivation_pt])/(totalbins - reactivation_pt))*reactivation_pt # y - mx, using reactivation_pt as the point from which y and x values are used
            rates["RapidExpo_2"].append((eachside_rate-rates["RapidExpo_2"][reactivation_pt])/(totalbins - reactivation_pt) * float(k) + intercept_RapidExpo_2)
            intercept_RapidExpo_10 = rates["RapidExpo_10"][reactivation_pt] - ((eachside_rate-rates["RapidExpo_10"][reactivation_pt])/(totalbins - reactivation_pt))*reactivation_pt # y - mx, using reactivation_pt as the point from which y and x values are used
            rates["RapidExpo_10"].append((eachside_rate-rates["RapidExpo_10"][reactivation_pt])/(totalbins - reactivation_pt) * float(k) + intercept_RapidExpo_10)
        #print i, rates["linear"]
        #print "xvals", x
        plt.figure()
        plt.plot(x,rates["Constant"],'r-', label="Constant")
        plt.plot(x,rates["Linear"],'g-', label="Linear")
        plt.plot(x,rates["Expo"],'b-', label="Exponential")
        #plt.plot(x,rates["RapidExpo_2"],'k-', label="Rapid Exponential-1")
        plt.plot(x,rates["RapidExpo_10"],'y-', label="Rapid Exponential")
        #plt.legend(loc=1, bbox_to_anchor=(1.0, 0.95), fontsize=9)
        xticktext = ["", "", "$\leq$0.5", "", "0.5-1", "", "1-1.5", "", "1.5-2", "", "2-2.5", "", "2.5-3", "", "3-3.5", "", "3.5-4", "", "4-4.5", "", "4.5-5", "", "5-5.5"]
        plt.xticks(range(0,origbins*2+1), xticktext, rotation=30) 
        plt.xlabel("Years in latency")
        plt.ylabel("New SNPs in year")
        #plt.ylim([0,24])
        plt.xlim([0,22])
        plt.tight_layout()
        #plt.title("Hypotheses of exponential decay of mutation incidence, $exp(\epsilon)$")
        plt.savefig("expgrowth_both_"+str(i)+"_6mo_bin.pdf")
        plt.close()
        binwise_rates_both.append(rates)

    #print(y1, "\n", y2, "\n", y4, "\n", y6)
    #return [y1, y2, y4, y5, y6]
    return binwise_rates_both
### end of getPoissonGrowth ##########

### start of getPoissonRateChange ##########
def getPoissonRateChange(rate): 
    x = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    y6 = []
    y7 = []
    factor = -1
    plt.figure()
    for i in range(0, len(N_monthwise)):
        x.append(i)
        y1.append(rate * np.exp(0))
        y2.append(-1 * rate/len(N_monthwise) * float(i) + rate)
        #y2.append(rate * np.exp(factor * float(i)/10.0))
        #y3.append(rate * np.exp(-1 * float(i)/3.0))
        y4.append(rate * np.exp(factor * i))
        y5.append(rate * np.exp(factor * i * 2))
        y6.append(rate * np.exp(factor * i * 10))
    #plt.plot(x,y1,'r-', label="$\epsilon$ = 0")
    #plt.plot(x,y2,'g-', label="$\epsilon$ = -x/10")
    #plt.plot(x,y4,'b-', label="$\epsilon$ = -x")
    #plt.plot(x,y5,'k-', label="$\epsilon$ = -x*2")
    #plt.plot(x,y6,'y-', label="$\epsilon$ = -x*10")
    plt.plot(x,y1,'r-', label="Constant")
    plt.plot(x,y2,'g-', label="Linear")
    plt.plot(x,y4,'b-', label="Exponential")
    #plt.plot(x,y5,'k-', label="Rapid Exponential-1")
    plt.plot(x,y6,'y-', label="Rapid Exponential")
    plt.legend(loc=1, bbox_to_anchor=(1.0, 0.95), fontsize=9)
    xticktext = ["$\leq$0.5","0.5-1", "1-1.5", "1.5-2", "2-2.5", "2.5-3", "3-3.5", "3.5-4", "4-4.5", "4.5-5", "5-5.5"]
    plt.xticks(range(0,len(N_monthwise)), xticktext, rotation=30) 
    plt.xlabel("Years in latency")
    plt.ylabel("New SNPs in year")
    plt.tight_layout()
    #plt.title("Hypotheses of exponential decay of mutation incidence, $exp(\epsilon)$")
    plt.savefig("expdecay.pdf")
    plt.close()

    #print(y1, "\n", y2, "\n", y4, "\n", y6)
    return [y1, y2, y4, y5, y6]
### end of getPoissonRateChange ##########

### start of murate_mean_ci ##########
def murate_mean_ci(yeargroup_mutperyear):
    # generate the figure that had murate wrt gen time
    yearkey = [[0,0,0,0,0,0,0,1,1,1,1],[2,2,2,3,3,3,3,3,3],[4,4,4,4,5]]
    gen = range(18,320,6)
    N = 0.973 * 4411532# 97.3% coverage for 4.41 MB genome of MTB. This is total number of sites in our analysis
    timeinhours = {} # key: year, value: that much time in hours. Year 0 => 0.5 years => 0.5 x 365 x 24 = hours in 0.5 years
    for t in range(0,6):
        timeinhours[t] = float((t+0.5)) * 365 * 24
    f, axarr = plt.subplots(1, 3, sharey=True, figsize=(8,6))
    axarr[0].set_ylabel("Mutation rate\n(mutations/bp/generation)")
    for i in range(len(yeargroup_mutperyear)):
        murates = {g : [] for g in gen} 
        for j in range(len(yeargroup_mutperyear[i])):
            year = yearkey[i][j]
            nsnps = yeargroup_mutperyear[i][j]
            t = timeinhours[year]
            print("nsnp:", nsnps, " year:", year)
            for g in gen:
                m = nsnps/(N*(t/g))
                murates[g].append(m) # append the mutation rate for this generation time

        yg_mean = []
        yg_95ci_upper = []
        yg_95ci_lower = []
        for g in gen:
            a = np.array(murates[g])
            a_mean = np.mean(a)
            output = [0.0,0.0]
            if a_mean != 0.0: # yeargrps 2-4 and 4-5 have 0 values, giving math error with st.norm.interval
                output = st.norm.interval(0.95, loc=a_mean, scale=np.std(a)/np.sqrt(len(a)))
            yg_mean.append(a_mean)
            yg_95ci_upper.append(output[1])
            yg_95ci_lower.append(output[0])
            if g == 18:
                print(yeargroup_mutperyear[i], a_mean, output)
        axarr[i].set_yscale('log')
        axarr[i].set_xlim([18,320])
        axarr[i].xaxis.set_ticks(np.arange(18, 320, 60))
        #axarr[i].axhline(y=mu_year01, c='gray', linestyle='--')
        axarr[i].plot(gen, yg_mean, lw = 1, color = '#539caf', alpha = 1)
        axarr[i].fill_between(gen, yg_95ci_lower, yg_95ci_upper, color = '#539caf', alpha = 0.4)
        #axarr[i].set_title("Years in latency: "+str(yg[0])+"-"+str(yg[1]+1), fontsize=10)
        break
    f.text(0.5, 0.04, 'Generation time (hours)', ha='center')
    plt.savefig("latency_yeargroups_yearwise_murates.pdf")
    exit()
                
### end of murate_mean_ci ##########

### start of drawFromPoisson ############
def drawFromPoisson(mrate, N_yearwise_orig, model, suffix, counter): # yearwise rate (= mean = expectation) of poission and num samples to be drawn
    # To increase sample size, use a multiplier for each year
    multiplier = 1
    N_yearwise_mod = [x * multiplier for x in N_yearwise_orig]
    n = sum(N_yearwise_mod)
    nleft = [0] + N_yearwise_mod # number of samples that are already drawn for the previous year
    yearwise_nmut = [] # list of lists: each row is the n draws from poission for that year using the rate for that year
    yearwise_avgmut = [] # avg num of mutations that year as per poisson distribution
    cumulative_nmut = [[0 for i in range(n)]] # cumulative number of mutations till that year
    cumulative_avgmut = [] # average num of mutations cumulative till that year
    ndraw = n
    indep_samples = [] # each row is an independent sample of mut accumulated uptil that year as per poisson dist with the given rate
    yearwise_mut_2_4_6 = [] # 1st row: mutations accumulated till year 2, 2nd row: mutations accumulated from year 2 to 4, 3rd row: mutations accumulated from year 4 to 6 years.
    indices_2_4_6 = [3, 7, len(mrate)-1] # the indices in mrate that correspond to 2 years, 4 years and the last time point
    cumulative_nmut_2_4_6 = [[0 for i in range(n)]] # cumulative number of mutations till years 2/4/last time point
    indep_samples_2_4_6 = []
    #print(mrate)
    for i in range(len(mrate)): # rate for that year
        m = mrate[i]
        ndraw = ndraw-nleft[i]
        mut = np.random.poisson(m, ndraw)
        yearwise_nmut.append(mut[:])
        yearwise_avgmut.append(np.mean(mut))
        cumulative = [x + y for x, y in zip(mut, cumulative_nmut[-1][0:ndraw])]
        indep_samples.append(cumulative[-nleft[i+1]:])
        cumulative_nmut.append(cumulative[:])
        #print("mut:", mut)
        #print("len_cumu_thisyear:", len(cumulative), "\tlen_mut:", len(mut), "\nndraw:", ndraw)
        #print("cumulative:", cumulative, "\nIndep:", indep_samples[-1])
        cumulative_2_4_6 = [x + y for x, y in zip(mut, cumulative_nmut_2_4_6[-1][0:ndraw])]
        cumulative_nmut_2_4_6.append(cumulative_2_4_6[:])
        indep_samples_2_4_6.extend(cumulative_2_4_6[-nleft[i+1]:])
        if i in indices_2_4_6: 
            #print("indep samples 2/4/6 at index ", i, "\n", indep_samples_2_4_6)
            temp = indep_samples_2_4_6[:]
            yearwise_mut_2_4_6.append(temp)
            indep_samples_2_4_6 = []
            cumulative_nmut_2_4_6 = [[0 for i in range(ndraw)]] # reset the mutations accumulated at the end of year 2 and 4
        cumulative_avgmut.append(np.mean(cumulative))
    #print("Final sets:", yearwise_mut_2_4_6)
    #murate_mean_ci(yearwise_mut_2_4_6)
    #exit()
    #print("Yearwise mutations: \n", yearwise_nmut)
    #print("\nCumulative mutations:\n", cumulative_nmut)
    #print("\nYearwise avg mut:", yearwise_avgmut)
    #print("\nCumulative avg mut:", cumulative_avgmut)
    fout = open("Data/"+model+"/"+model+"_poisson_"+str(n)+"samples_"+str(counter)+suffix+".txt", 'w')
    for i in range(len(indep_samples)):
        for j in range(len(indep_samples[i])):
            fout.write(str(i)+","+str(indep_samples[i][j])+"\n")
    fout.close()
    #exit()
    plt.plot(range(len(mrate)), mrate, label="Rate")
    plt.plot(range(len(mrate)), yearwise_avgmut, label="Avg SNPs/year")
    plt.plot(range(len(mrate)), cumulative_avgmut, label="Avg. Cumulative SNPs")
    plt.xlabel("Years in latency")
    plt.ylabel("Num. SNPs")
    plt.legend(loc=2, ncol=3)
    plt.savefig("Figures/"+model+"/"+model+"_mut_trends_"+str(n)+"samples_"+str(counter)+".pdf")
    plt.close()
### end of drawFromPoisson ############

### start of plotPvals #############
def plotPvals(fname, modelcolor): # read p-values of stat test that compare true data to sampled data acc to diff models
    fh = open(fname, 'r')
    plt.figure()
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split(",")
        modelname = contents[0]
        modelpvals = [float(v) for v in contents[1:]]
        histvals, binedges = np.histogram(np.log10(modelpvals), bins=10) # take log10 vals of p-values
        print(modelname, "\n", len(histvals), "\n", len(binedges[:-1]))
        plt.plot(binedges[:-1], histvals, modelcolor[modelname], label=modelname)
    plt.xlabel("$log_{10}(p-value)$")
    plt.ylabel("Frequency")
    plt.legend(loc=2, ncol=3)
    plt.savefig("latencymodels_pvalhist.pdf")
    plt.close()
    fh.close()
### end of plotPvals #############

### start of plotBoxplot #############
def plotBoxplot(fname, modelxtext, prefix): # read p-values of stat test that compare true data to sampled data acc to diff models
    fh = open(fname, 'r')
    boxplotdata = []
    xtext = []
    plt.figure()
    print("\nMedian p-values:\n")
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split(",")
        modelname = contents[0]
        if modelname.startswith("RapidExpo_2"): # keeping only one rapid-exponential model, RapidExpo_10
            continue # ignore RapidExpo_2 model data
        modelpvals = [np.log10(float(v)) for v in contents[1:]]
        print(modelname, np.median(modelpvals))
        boxplotdata.append(modelpvals)
        xtext.append(modelxtext[modelname])
    plt.boxplot(boxplotdata, 0, '+')
    plt.ylabel("$log_{10}(p-value)$")
    plt.xticks(range(1,len(xtext)+1),xtext)
    plt.xticks(fontsize=9)
    plt.axhline(y=-2, color='gray', linestyle='--')
    plt.axhline(y=-1.3, color='black', linestyle='--')
    plt.savefig(prefix+"_latencymodels_pvalboxplot.pdf")
    plt.close()
    fh.close()
### end of plotBoxplot #############

### start of both_poisson ##########
def both_poisson(latencyrates, N, nrep, modelname):
    indexcol = []
    counter = 0
    for nval in N:
        indexcol.extend([counter] * nval)
        counter += 1
    for i in range(nrep): # for each replicate
        allcumul = {"Constant":[], "Linear":[], "Expo":[], "RapidExpo_2":[], "RapidExpo_10":[]}
        for j in range(len(latencyrates)): # each bin (6-month window) has its own rates
            for m in modelname:
                rate = latencyrates[j][m]
                #n_thisrate = [0] * (len(rate)-1)
                #n_thisrate.append(N[j])
                ndraw = N[j] # number of latency pairs in this bin
                #print i, "\tbin:", j, m, rate, ndraw
                #drawFromPoisson(rate, n_thisrate, m+"_both", "_bin"+str(j), i) # early latency model only
                mut = []
                cumulative = [0] * ndraw

                for k in range(len(rate)): # rate for that year
                    mut = np.random.poisson(rate[k], ndraw)
                    cumulative = [x + y for x, y in zip(mut, cumulative)]
                    #print(k, "mut:", mut)
                    #print("cumulative:", cumulative)
                allcumul[m].extend(cumulative)
        for m in allcumul:
            model = m+"_both"
            fout = open("Data/"+model+"/"+model+"_poisson_"+str(len(allcumul[m]))+"samples_"+str(i)+".txt", 'w')
            for p in range(len(allcumul[m])): 
                fout.write(str(indexcol[p])+","+str(allcumul[m][p])+"\n")
            fout.close()
        #exit()
### end of both_poisson ##########

### MAIN ###
modelname = ["Constant", "Linear", "Expo", "RapidExpo_2", "RapidExpo_10"]
modelcolor = {"Constant":'r-', "Linear":'g-', "Expo":'b-', "RapidExpo_2":'k-', "RapidExpo_10":'y-', "RapidExpo_100":'c-'}
modelxtext = {"Constant":"Constant", "Linear":"Linear", "Expo":"Exponential", "RapidExpo_2":"Rapid\nExponential-1", "RapidExpo_10":"Rapid\nExponential"}
modelxtext_reactivation = {"Constant_reactivation":"Constant", "Linear_reactivation":"Linear", "Expo_reactivation":"Exponential", "RapidExpo_2_reactivation":"Rapid\nExponential-1", "RapidExpo_10_reactivation":"Rapid\nExponential"}
modelxtext_both = {"Constant_both":"Constant", "Linear_both":"Linear", "Expo_both":"Exponential", "RapidExpo_2_both":"Rapid\nExponential-1", "RapidExpo_10_both":"Rapid\nExponential"}
#latencyrates = getPoissonReactivation(rate)
#latencyrates_both = getPoissonGrowthBoth(avgsnp_monthwise)
#latencyrates = getPoissonRateChange(rate)
#exit()
#print len(latencyrates)
#plotPvals("pvalues_ks2s2d_latencymodels.txt", modelcolor)
plotBoxplot("pvalues_ks2s2d_monthwise_latencymodels.txt", modelxtext, "early_only")
plotBoxplot("pvalues_ks2s2d_monthwise_latencymodels_reactivation.txt", modelxtext_reactivation, "reactivation")
#plotBoxplot("pvalues_ks2s2d_monthwise_latencymodels_both.txt", modelxtext_both, "both")
exit()
Nrand = 1000 # num of independent samples to be drawn
#both_poisson(latencyrates_both, N_monthwise, Nrand, modelname) # model both early decline and reactivation
#exit()
for i in range(len(latencyrates)):
    rate = latencyrates[i]
    print("Model: ", modelname[i], "\tYearwise rates:", rate)
#    if modelname[i] != "RapidExpo_10":
#        continue
    for k in range(Nrand):
        if k % 100 == 0:
            print k
        drawFromPoisson(rate, N_monthwise, modelname[i], "", k) # early latency model only
        #drawFromPoisson(rate, N_monthwise, modelname[i]+"_reactivation", "", k) # reactivation model only
