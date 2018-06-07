#!/usr/bin/python

import scipy.stats as ss

early_age = [24, 39, 18, 19, 17, 36, 19, 42, 19, 6, 18]
late_age = [27, 72, 15, 16, 16, 42, 40, 24, 48, 19, 21, 23, 32, 21]
statval_age, pval_age = ss.ranksums(early_age, late_age)
print ("Pval for age comparison: ", pval_age)

bmi_early = [20.96, 17.21, 17.56, 23.81, 19.08, 20.94]
bmi_late = [18.8, 19.3, 27.3, 23.53, 17.5, 23.97, 24.3, 21.3, 19.37, 17.92, 18.21]
statval_bmi, pval_bmi = ss.ranksums(bmi_early, bmi_late)
print ("Pval for bmi comparison: ", pval_bmi)

ppdmax_early = [0, 13, 19, 11, 20, 16, 20, 13, 8, 9]
ppdmax_late = [15, 11, 9, 15, 21, 17, 15, 18]
statval_ppdmax, pval_ppdmax = ss.ranksums(ppdmax_early, ppdmax_late)
print ("Pval for ppdmax comparison: ", pval_ppdmax)


early_gender = [7, 4]
late_gender = [5, 9]
oddsratio_gender, pval_gender = ss.fisher_exact([early_gender, late_gender]) # requires 2x2 matrix
print ("Pval for gender comparison: ", pval_gender)

early_bcgscar = [10, 1]
late_bcgscar = [11, 2]
oddsratio_bcgscar, pval_bcgscar = ss.fisher_exact([early_bcgscar, late_bcgscar])
print ("Pval for bcgscar comparison: ", pval_bcgscar)

early_hiv = [0, 6]
late_hiv = [0, 12]
oddsratio_hiv, pval_hiv = ss.fisher_exact([early_hiv, late_hiv])
print ("Pval for hiv comparison: ", pval_hiv)

early_alcohol = [1, 4] # same numbers are for smoking as well
late_alcohol = [1, 5]
oddsratio_alcohol, pval_alcohol = ss.fisher_exact([early_alcohol, late_alcohol])
print ("Pval for alcohol and smoking comparison: ", pval_alcohol)

early_comorb = [0, 5]
late_comorb = [2, 5]
oddsratio_comorb, pval_comorb = ss.fisher_exact([early_comorb, late_comorb])
print ("Pval for comorb comparison: ", pval_comorb)
