# NCFig3BFig4.R
# generate Figure 4 and Figure 3B

rm(list=ls())
#install.packages("sandwich")
#install.packages("aod")
#update.packages()

require(sandwich)
require(aod)

#-----------------------------------------------------------------------
# Read in the data and generate variables
snp.data=read.csv(file="C:/mydoc/Rutgers/Aditi/latency_snplist.csv")
#remove the summary stats that were at the bottom of the file
snp.data=snp.data[1:24,]
#original sort order
origorder=1:24

# create groupings of year of diagnosis of secondary case
yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,72), include.lowest=TRUE, right=T, 
         label=c("[0,2]","(2,6]"))

# add new breakdown of year of diagnosis to the dataframe
snp.data=data.frame(snp.data,yrgp,origorder)
# put the data in order of months
ord=order(snp.data$Months.between.diagnosis)
snp.data=snp.data[ord,]
rm(ord,yrgp,origorder)

#--------------------------------------------------------------------------------------
# range of possible generation times
gentime=seq(1,340) # from 18 to 318 hours (expand by +22/-17 so it reaches the edges of plot)
# MTB genome details and coverage
gensize=4411532 # assumed length of genome
coverage=0.973 #97.3% coverage of the genome
# factor to get the rate per bp*hour (assume 30.4 days/month*24 hours/day) (note I did not take the log of this quantity)
fctr =coverage*gensize*30.4*24/gentime
# print out the fctr when the generation time is 18 hours (for checking)
fctr[gentime==18]
#[1] 173986116

#--------------------------------------------------------------------------------------
#fit the Poisson model with Months btwn diagnosis as a linear term (note: no offset used)
M1=glm(formula=X..SNPs.different~Months.between.diagnosis , family="poisson", data=snp.data)

#reduced Poisson model (intercept only)
M0=glm(formula=X..SNPs.different~1 , family="poisson", data=snp.data)

#test whether Month.between.diagnosis is significant using LRT
anova(M1,M0)
pchisq(2*(logLik(M1) - logLik(M0)), df = 1, lower.tail = FALSE)
#'log Lik.' 0.8821014 (df=2) -> p=0.88
#
#robust variance/SE
cov.M1=vcovHC(M1,type="HC0")
std.err.M1=sqrt(diag(cov.M1))

# wald test using robust variance
l=matrix(c(0,1),nrow=1, ncol=2)
wt=wald.test(b=coef(M1), Sigma=cov.M1, L=l)
wt
#Wald test:
#----------
#  
#  Chi-squared test:
#  X2 = 0.017, df = 1, P(> X2) = 0.9

# this is the poisson estimate and 95% CI overall obtained from the reduced (intercept only) model
est.M0=c(Estimate=coef(M0),confint(M0))
# calculate robust variance/SE and estimate with 95% CIs using them 
cov.M0=vcovHC(M0,type="HC0")
std.err.M0=sqrt(cov.M0)

r.est.M0 <- c(Estimate= coef(M0), "Robust SE" = std.err.M0,
               "Pr(>|z|)" = 2 * pnorm(abs(coef(M0)/std.err.M0), lower.tail=FALSE),
               LL = coef(M0) - 1.96 * std.err.M0,
               UL = coef(M0) + 1.96 * std.err.M0)
r.est.M0
#Estimate.(Intercept)            Robust SE             Pr(>|z|)                   LL                   UL 
#         0.810930216          0.281152342          0.003922737          0.259871627          1.361988806 

# Poisson estimate of SNPs with robust 95% CIs
poissrate<- c(Estimate= exp(coef(M0)), 
              LL = exp(coef(M0) - 1.96 * std.err.M0),
              UL = exp(coef(M0) + 1.96 * std.err.M0))
poissrate
#Estimate.(Intercept)       LL.(Intercept)       UL.(Intercept) 
#            2.250000             1.296764             3.903950
#-----------------------------------------------------------------------------
# calculate robust SE and Poisson estimate use robust SE to calculate 95% CIs from them M1 model with months as a predictor (see above)
cov.M1=vcovHC(M1,type="HC0")
std.err.M1=sqrt(diag(cov.M1))

r.est.M1 <- cbind(Estimate= coef(M1), "Robust SE" = std.err.M1,
               "Pr(>|z|)" = 2 * pnorm(abs(coef(M1)/std.err.M1), lower.tail=FALSE),
               LL = coef(M1) - 1.96 * std.err.M1,
               UL = coef(M1) + 1.96 * std.err.M1)
r.est.M1
#                             Estimate   Robust SE    Pr(>|z|)          LL         UL
#(Intercept)               0.839811612 0.322916368 0.009303124  0.20689553 1.47272769
#Months.between.diagnosis -0.001023528 0.007765929 0.895144660 -0.01624475 0.01419769

#note: the above shows p=0.895144660; therefore we reject the hypothesis that SNPs differ by Months.between.diagnosis

#----------------------------------------------------------------
# Estimate a rate per bp*generation
# if values of Months between diagnosis = 0 set to 1 since log(0) is undefined

# make a new dataframe for this analysis so we are not changing the data in the original dataframe
snp.data1=snp.data
# if months btwn diagnoses for the TB pair is 0, then set to it 1 since you cannot take the log of 0
snp.data1$Months.between.diagnosis[snp.data1$Months.between.diagnosis==0]=1
snp.data1$Months.between.diagnosis

#generate logtime 
logtime=log(snp.data1$Months.between.diagnosis)

#generate log number of generations*bp assuming generation time is 18 hours
lognbgen=log(snp.data1$Months.between.diagnosis*coverage*gensize*30.4*24/18)

#add new variables to dataframe
snp.data1=data.frame(snp.data1,logtime,lognbgen)
rm(logtime)
rm(lognbgen)
#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis was [0,2] years using Months between
# diagnoses as the time at risk (by using the log of the time as an offset)
M0to2yr=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis),
              family="poisson", data=snp.data1, subset=yrgp=="[0,2]")

# this is the poisson estimate and model based 95% CI 
summary(M0to2yr)
est.M0to2yr=c(Estimate=coef(M0to2yr),confint(M0to2yr))

# calculate robust standard errors
cov.M0to2yr=vcovHC(M0to2yr,type="HC0")
std.err.M0to2yr=sqrt(cov.M0to2yr)

r.est.M0to2yr <- c(Estimate= coef(M0to2yr), "Robust SE" = std.err.M0to2yr,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(M0to2yr)/std.err.M0to2yr), lower.tail=FALSE),
                   LL = coef(M0to2yr) - 1.96 * std.err.M0to2yr,
                   UL = coef(M0to2yr) + 1.96 * std.err.M0to2yr)
r.est.M0to2yr
#Estimate.(Intercept)            Robust SE             Pr(>|z|)                   LL                   UL 
#-1.927892e+00         4.143756e-01         3.278993e-06        -2.740068e+00        -1.115715e+00 

# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was 0-2 years
poissrate0to2yr<- c(Estimate= exp(coef(M0to2yr)), 
                    LL = exp(coef(M0to2yr) - 1.96 * std.err.M0to2yr),
                    UL = exp(coef(M0to2yr) + 1.96 * std.err.M0to2yr))
poissrate0to2yr
#Estimate.(Intercept)                   LL                   UL 
#          0.14545455           0.06456597           0.32768074 

# check to see if my estimate from glm is correct by calc SNP/months
sum(snp.data1$X..SNPs.different[snp.data$yrgp=="[0,2]"])
#16
sum(snp.data1$Months.between.diagnosis[snp.data$yrgp=="[0,2]"])
#110 months
##estimate
#16/110
#[1] 0.1454545
# matches above estimate from glm

# now use the range of generation times to obtain estimates of rates and robust SE based 95% CIs
rate.M0to2yr=poissrate0to2yr[1]/fctr
rate.M0to2yr.LL=poissrate0to2yr[2]/fctr
rate.M0to2yr.UL=poissrate0to2yr[3]/fctr

# rate and 95% CI for gentime=18
cbind(rate.M0to2yr,rate.M0to2yr.LL,rate.M0to2yr.UL)[gentime==18,]
#rate.M0to2yr rate.M0to2yr.LL rate.M0to2yr.UL 
#8.360124e-10    3.710984e-10    1.883373e-09 

#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis of 2+ years

M2plusyr=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis),
            family="poisson", data=snp.data1, subset=yrgp=="(2,6]")
# this is the poisson estimate and model based 95% CI 
summary(M2plusyr)
est.M2plusyr=c(Estimate=coef(M2plusyr),confint(M2plusyr))
est.M2plusyr
# let's now do robust standard errors
cov.M2plusyr=vcovHC(M2plusyr,type="HC0")
std.err.M2plusyr=sqrt(cov.M2plusyr)

# unexponentiated point estimate and standard error
coef(M2plusyr)
#(Intercept) 
#-2.7133 
std.err.M2plusyr
#(Intercept)   0.4110211

r.est.M2plusyr <- c(Estimate= coef(M2plusyr), "Robust SE" = std.err.M2plusyr,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(M2plusyr)/std.err.M2plusyr), lower.tail=FALSE),
                   LL = coef(M2plusyr) - 1.96 * std.err.M2plusyr,
                   UL = coef(M2plusyr) + 1.96 * std.err.M2plusyr)
r.est.M2plusyr
# Estimate.(Intercept)            Robust SE             Pr(>|z|)                   LL                   UL 
#        -2.713300e+00         4.110211e-01         4.073934e-11        -3.518901e+00        -1.907698e+00 

# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was >4 years
poissrate2plusyr<- c(Estimate= exp(coef(M2plusyr)), 
                    LL = exp(coef(M2plusyr) - 1.96 * std.err.M2plusyr),
                    UL = exp(coef(M2plusyr) + 1.96 * std.err.M2plusyr))
poissrate2plusyr
#Estimate.(Intercept)                   LL                   UL 
#          0.06631763           0.02963199           0.14842163 

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M2plusyr=poissrate2plusyr[1]/fctr
rate.M2plusyr.LL=poissrate2plusyr[2]/fctr
rate.M2plusyr.UL=poissrate2plusyr[3]/fctr

# rate and 95% CI for gentime=18
cbind(rate.M2plusyr,rate.M2plusyr.LL,rate.M2plusyr.UL)[gentime==18,]
#rate.M2plusyr rate.M2plusyr.LL rate.M2plusyr.UL 
#3.811662e-10     1.703124e-10     8.530659e-10 

#=============================================================
# create a PDF file with the two plots side by side
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\NCFig4.pdf",height=6,width=10)
par(mfrow=c(1,2), mar=c(4,3,2,0)+0.1,oma=c(0,5,1,1))
# this is the yaxis limits for the graph
ylim=c(1e-11,1e-6)
#----------------------------------------------------------
# 0-2 yrs
plot(gentime, rate.M0to2yr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: [0, 2] Years",
     xlab="Generation time (hours)", 
     ylab="", 
     col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1,cex.lab=1.5)
#ylabel
title(ylab="Mutation Rate (bp/generation)",outer=T, line=0,cex.lab=1.5)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M0to2yr.LL,rev(rate.M0to2yr.UL)),
        col=adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)
# put y-axis tick marks w/o label
axis(side=2,at=c(seq(1:9)*1e-11,
                 seq(1:9)*1e-10,seq(1:9)*1e-9,seq(1:9)*1e-8,seq(1:9)*1e-7,seq(1:9)*1e-6),
     labels=F,tcl=-.2)
# add y-axis labels only at 10^super
super=c(-11,-10,-9,-8,-7,-6)
yat=10^super
ylab=parse(text=paste(rep("10^",length(super)),super,sep=""))
axis(side=2,at=yat, label=ylab, las=2,tcl=-.3, lwd=1,adj=.5) #tcl is used to shorten tick marks
# put x-axis tick marks and labels 
xat=18+60*(0:5)
axis(side=1,at=xat, label=xat, las=1,lwd=1, adj=.5)

#add horizontal line at rate for 0-2 year interval rate for generation time of 18 hrs
abline(h=rate.M0to2yr[gentime==18],lty=2,col="grey")
abline(v=18,lty=2,col="grey")
legend(168,1.5*1e-11,legend=expression(paste("rate=8.36×10"^"-10")),lty=2,bty="n",col="grey",cex=1)

#----------------------------------------------------------
plot(gentime, rate.M2plusyr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: (2, 6] Years",
     xlab="Generation time (hours)", 
    # ylab="Mutation Rate (mutations/bp/generation)", 
    ylab="" ,
    col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1,cex.lab=1.5)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M2plusyr.LL,rev(rate.M2plusyr.UL)),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)
# put y-axis tick marks w/o label
axis(side=2,at=c(seq(1:9)*1e-11,seq(1:9)*1e-10,seq(1:9)*1e-9,seq(1:9)*1e-8,seq(1:9)*1e-7,seq(1:9)*1e-6),
     labels=F,tcl=-.2)
# add y-axis labels only at 10^super
super=c(-11,-10,-9,-8,-7,-6)
yat=10^super
ylab=parse(text=paste(rep("10^",length(super)),super,sep=""))
axis(side=2,at=yat, label=ylab, las=2,tcl=-.3, lwd=1,adj=.5) #tcl is used to shorten tick marks

# put x-axis tick marks and labels
xat=18+60*(0:5)
axis(side=1,at=xat, label=xat, las=1,lwd=1, adj=.5)
#add horizontal line at rate for 0-2 year interval rate for generation time of 18 hrs
abline(h=rate.M0to2yr[gentime==18],lty=2,col="grey")
abline(v=18,lty=2,col="grey")
legend(168,1.5*1e-11,legend=expression(paste("rate=8.36×10"^"-10")),lty=2,bty="n",col="grey",cex=1)
#----------------------------------------------------------
dev.off()
#============================================================================================
# model with offset and 2 latency period categories
M2gp=glm(formula=X..SNPs.different~yrgp , offset= log(Months.between.diagnosis), family="poisson", data=snp.data1)
summary(M2gp)

cov.M2gp=vcovHC(M2gp,type="HC0")
std.err.M2gp=sqrt(diag(cov.M2gp))

#reduced Poisson model (intercept only)
M1gp=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis), family="poisson", data=snp.data)

#test whether yrgp is significant using LRT test
anova(M2gp,M1gp,test="Chisq")
as.numeric(pchisq(2*(logLik(M2gp) - logLik(M1gp)), df = 1, lower.tail = FALSE))
#[1] 0.01315706
# p=0.01315706

# test whether yrgp is significant using the robust SE in Wald test (more conservative)
l=matrix(c(0,1),nrow=1, ncol=2)
wt=wald.test(b=coef(M2gp), Sigma=cov.M2gp, L=l)
wt
#Wald test:
#  ----------
#  
#  Chi-squared test:
#  X2 = 1.8, df = 1, P(> X2) = 0.18

#================================================================================
# For Fig 3B
# Fit model SNP rate per bp*generation with months betwen IC-HHC diagnosis as a continuous variable for 18 hour generation time
fctr18=coverage*gensize*30.4*24/18
M10=glm(formula=X..SNPs.different~Months.between.diagnosis,family="poisson",
        offset=log(Months.between.diagnosis*fctr18),data=snp.data1)

# calculate robust standard errors
cov.M10=vcovHC(M10,type="HC0")
std.err.M10=sqrt(diag(cov.M10))

r.est.M10 <- cbind(Estimate= coef(M10), "Robust SE" = std.err.M10,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(M10)/std.err.M10), lower.tail=FALSE),
                   LL = coef(M10) - 1.96 * std.err.M10,
                   UL = coef(M10) + 1.96 * std.err.M10)
r.est.M10
#                                      Robust SE     Pr(>|z|)           LL           UL
#(Intercept)              -19.64359757 0.44912212 0.000000e+00 -20.52387693 -18.76331822
#Months.between.diagnosis  -0.05320629 0.01247237 1.990696e-05  -0.07765213  -0.02876044

#robust variance
cov.M10=vcovHC(M10,type="HC0")
std.err.M10=sqrt(diag(cov.M10))

# wald test using robust variance
l=matrix(c(0,1),nrow=1, ncol=2)
wt=wald.test(b=coef(M10), Sigma=cov.M10, L=l)
wt
#Wald test:
#  ----------
#
#  Chi-squared test:
#  X2 = 18.2, df = 1, P(> X2) = 2e-05
# p=2e-05 which is <0.001

pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\NCFig3B.pdf", height=7, width=7)
#par (mar=c(4,6,1,1)+0.1)
par(mar=c(6,6,5,2)+0.1)
# in the following plot the IC-HHC pair with 0 months between diagnoses is plotted
#   at month 0 but the vertical axis is plotted assuming that it occurred at month 1 (this is
#   because it should be divided by the months which is undefined

plot(snp.data$Months.between.diagnosis, 
        snp.data$X..SNPs.different/snp.data$Months.between.diagnosis/fctr18,
     xlim=c(1,64),xlab="Months between Index Case and HHC TB Diagnosis", 
     ylim=c(0,2.5e-8),ylab="Mutation Rate (mutations/bp/generation)\n\n", 
     cex.lab=1, yaxt="n", xaxt="n" , pch=20)
     #cex=.8 )
points(snp.data$Months.between.diagnosis[1], 
         snp.data1$X..SNPs.different[1]/snp.data1$Months.between.diagnosis[1]/fctr18,pch="^",
       cex=0.8)

# create a dataset of 0.1 month increments so we can plot a smooth line
newdata=data.frame(Months.between.diagnosis=seq(1,63,by=.1))
pred=predict(M10,newdata=newdata,type="response")/(newdata$Months.between.diagnosis*fctr18)
lines(newdata$Months.between.diagnosis,pred)
# add y-axis labels (I put the value of 0 separately so it is 0 not 0e0)
attk=seq(1e-9,1e-7,by=1e-9)
axis(side=2,at=attk,labels=F,tcl=-.2,las=2,cex.axis=1)
attk=seq(5e-9,1e-7,by=5e-9)
axis(side=2,at=attk,labels=T,tcl=-.5,las=2,cex.axis=1)
axis(side=2,at=0,labels=T,tcl=-.2,las=2,cex.axis=1)
# add x-axis label
axis(side=1,at=seq(0,66,by=6),labels=T,tcl=-.5,cex.axis=1)
#axis(side=1,at=seq(0,65,by=5),labels=F,tcl=-.2,cex.axis=.7)
#axis(side=1,at=seq(0,65,by=10),labels=T,tcl=-.5,cex.axis=.7)

dev.off()
#-------------------------------------------------------------------------------------
#output observed and predicted rates
out=cbind(snp.data$Index.case, snp.data$HHC ,snp.data$Months.between.diagnosis,snp.data$X..SNPs.different,snp.data1$X..SNPs.different/snp.data1$Months.between.diagnosis/fctr18,
      predict(M10,snp.data1,type="response")/(snp.data1$Months.between.diagnosis*fctr18) )
out=out[order(snp.data$origorder),]

dimnames(out)[[2]]=c("Index.case","HHC","Months latency","SNPs","Observed Rate","Predicted Rate")
out

write.table(out,sep=",",file="ObsPredRates.csv",row.names=FALSE)

