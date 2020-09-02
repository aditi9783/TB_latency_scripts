# NCSuppFig8.R 
# Generate Supplementary Figure 8

rm(list=ls())
#install.packages("sandwich")
#install.packages("aod")
#update.packages()
require(sandwich)
require(aod)

#
#-----------------------------------------------------------------------
# Read in the data and generate variables
snp.data=read.csv(file="C:/mydoc/Rutgers/Aditi/latency_snplist.csv")
#remove the summary stats that were at the bottom of the file
snp.data=snp.data[1:24,]
# put the data in order of months -- useful later when plotting
ord=order(snp.data$Months.between.diagnosis)
snp.data=snp.data[ord,]

# create a coarser breakdown of year of diagnosis of secondary case and add to data frame
yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,48,72), include.lowest=TRUE, right=T, 
           label=c("[0,2]","(2,4]","(4,6]"))

# add new breakdown of year of diagnosis to the dataframe
snp.data=data.frame(snp.data,yrgp)
rm(yrgp)
names(snp.data)
summary(snp.data)
#--------------------------------------------------------------------------------------
# set parameters that will be used later

# range of possible generation times
gentime=seq(1,340) # from 18 to 318 hours (+22/-17 so it reaches the edges of plot)
# genome details and coverage
gensize=4411532 # assumed length of genome
coverage=0.973 #97.3% coverage of the genome
# factor to get the rate per bp*hour (assume 30.4 days/month*24 hours/day) (note I did not take the log of this quantity)
fctr =coverage*gensize*30.4*24/gentime
# print out the fctr when the generation time is 18 hours (for checking)
fctr[gentime==18]
#[1] 173986116
#----------------------------------------------------------------
#  make a new dataframe for this analysis so we are not changing the data in the original dataframe
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
#fit the Poisson model for HHCs whose time of diagnosis was [0,2] years (or [0,24] months)
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

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M0to2yr=poissrate0to2yr[1]/fctr
rate.M0to2yr.LL=poissrate0to2yr[2]/fctr
rate.M0to2yr.UL=poissrate0to2yr[3]/fctr

# rate and 95% CI for gentime=18
cbind(rate.M0to2yr,rate.M0to2yr.LL,rate.M0to2yr.UL)[gentime==18,]
#rate.M0to2yr rate.M0to2yr.LL rate.M0to2yr.UL 
#8.360124e-10    3.710984e-10    1.883373e-09 

#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis of (2,4] years (or (24,48] months)

M2to4yr=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis),
            family="poisson", data=snp.data1, subset=yrgp=="(2,4]")

# this is the poisson estimate and model based 95% CI 
summary(M2to4yr)
est.M2to4yr=c(Estimate=coef(M2to4yr),confint(M2to4yr))

# calculate robust standard errors
cov.M2to4yr=vcovHC(M2to4yr,type="HC0")
std.err.M2to4yr=sqrt(cov.M2to4yr)

r.est.M2to4yr <- c(Estimate= coef(M2to4yr), "Robust SE" = std.err.M2to4yr,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(M2to4yr)/std.err.M2to4yr), lower.tail=FALSE),
                   LL = coef(M2to4yr) - 1.96 * std.err.M2to4yr,
                   UL = coef(M2to4yr) + 1.96 * std.err.M2to4yr)
r.est.M2to4yr 
#Estimate.(Intercept)            Robust SE             Pr(>|z|)                   LL                   UL 
#-2.326538e+00         4.401259e-01         1.249691e-07        -3.189185e+00        -1.463892e+00 


# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was >2-4 years
poissrate2to4yr<- c(Estimate= exp(coef(M2to4yr)), 
                    LL = exp(coef(M2to4yr) - 1.96 * std.err.M2to4yr),
                    UL = exp(coef(M2to4yr) + 1.96 * std.err.M2to4yr))
poissrate2to4yr
#Estimate.(Intercept)                   LL                   UL 
#          0.09763314           0.04120543           0.23133428 

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M2to4yr=poissrate2to4yr[1]/fctr
rate.M2to4yr.LL=poissrate2to4yr[2]/fctr
rate.M2to4yr.UL=poissrate2to4yr[3]/fctr

# rate and 95% CI for gentime=18
cbind(rate.M2to4yr,rate.M2to4yr.LL,rate.M2to4yr.UL)[gentime==18,]
#rate.M2to4yr rate.M2to4yr.LL rate.M2to4yr.UL 
#5.611548e-10    2.368317e-10    1.329613e-09 

#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis of 4+ years (or (48+ months)

M4plusyr=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis),
            family="poisson", data=snp.data1, subset=yrgp=="(4,6]")
# this is the poisson estimate and model based 95% CI 
summary(M4plusyr)
est.M4plusyr=c(Estimate=coef(M4plusyr),confint(M4plusyr))
est.M4plusyr
# let's now do robust standard errors
cov.M4plusyr=vcovHC(M4plusyr,type="HC0")
std.err.M4plusyr=sqrt(cov.M4plusyr)

r.est.M4plusyr <- c(Estimate= coef(M4plusyr), "Robust SE" = std.err.M4plusyr,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(M4plusyr)/std.err.M4plusyr), lower.tail=FALSE),
                   LL = coef(M4plusyr) - 1.96 * std.err.M4plusyr,
                   UL = coef(M4plusyr) + 1.96 * std.err.M4plusyr)
r.est.M4plusyr
#Estimate.(Intercept)            Robust SE             Pr(>|z|)                   LL                   UL 
#-3.850148e+00         4.224141e-01         7.894061e-20        -4.678079e+00        -3.022216e+00 

# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was >4 years
poissrate4plusyr<- c(Estimate= exp(coef(M4plusyr)), 
                    LL = exp(coef(M4plusyr) - 1.96 * std.err.M4plusyr),
                    UL = exp(coef(M4plusyr) + 1.96 * std.err.M4plusyr))
poissrate4plusyr
#Estimate.(Intercept)                   LL                   UL 
#         0.021276596          0.009296854          0.048693197 

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M4plusyr=poissrate4plusyr[1]/fctr
rate.M4plusyr.LL=poissrate4plusyr[2]/fctr
rate.M4plusyr.UL=poissrate4plusyr[3]/fctr

# rate and 95% CI for gentime=18
cbind(rate.M4plusyr,rate.M4plusyr.LL,rate.M4plusyr.UL)[gentime==18,]
#rate.M4plusyr rate.M4plusyr.LL rate.M4plusyr.UL 
#1.222890e-10     5.343446e-11     2.798683e-10 

#=============================================================
# create a PDF file with all three plots side by side
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\NCSuppFig8.pdf",height=6,width=10)
par(mfrow=c(1,3), mar=c(4,3,2,0)+0.1,oma=c(0,5,1,1))
# this is the yaxis limits for the graph
ylim=c(1e-11,1e-6)
#----------------------------------------------------------
# 0-2 yrs
plot(gentime, rate.M0to2yr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: [0, 2] Years",
     xlab="Generation time (hours)", 
#     ylab="Mutation Rate (mutations/bp/generation)", 
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
#add vertical line
abline(v=18,lty=2,col="grey")

#----------------------------------------------------------
# 2-4 yrs
plot(gentime, rate.M2to4yr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: (2, 4] Years",
     xlab="Generation time (hours)", 
    # ylab="Mutation Rate (mutations/bp/generation)", 
     ylab="" ,
     col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1,cex.lab=1.5)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M2to4yr.LL,rev(rate.M2to4yr.UL)),
        col=adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)
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
#add vertical line
abline(v=18,lty=2,col="grey")

#----------------------------------------------------------
plot(gentime, rate.M4plusyr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: (4, 6] Years",
     xlab="Generation time (hours)", 
    # ylab="Mutation Rate (mutations/bp/generation)", 
    ylab="" ,
    col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1,cex.lab=1.5)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M4plusyr.LL,rev(rate.M4plusyr.UL)),
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
#add vertical line
abline(v=18,lty=2,col="grey")
legend(168,1.5*1e-11,legend=expression(paste("rate=8.36×10"^"-10")),lty=2,bty="n",col="grey",cex=1)
#----------------------------------------------------------
dev.off()
#============================================================================================
# model with offset and 3 latency period categories
M3gp=glm(formula=X..SNPs.different~yrgp , offset= log(Months.between.diagnosis), family="poisson", data=snp.data1)
summary(M3gp)

cov.M3gp=vcovHC(M3gp,type="HC0")
std.err.M3gp=sqrt(diag(cov.M3gp))

#reduced Poisson model (intercept only)
M1gp=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis), family="poisson", data=snp.data1)

#test whether yrgp is significant using LRT test
anova(M3gp,M1gp,test="Chisq")
as.numeric(pchisq(2*(logLik(M3gp) - logLik(M1gp)), df = 2, lower.tail = FALSE))
#[1] 3.896625e-05
# p=0.00003896625

# test whether yrgp is significant using the robust SE in Wald test (more conservative)
l=rbind(c(0,1,0),c(0,0,1))
wt=wald.test(b=coef(M3gp), Sigma=cov.M3gp, L=l)
wt
#Wald test:
#----------
#  
# Chi-squared test:
# X2 = 11.6, df = 2, P(> X2) = 0.003
# p=0.003

