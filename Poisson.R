# Poisson regression analyses
#install.packages("sandwich")
#install.packages("MASS")
#update.packages()
require(ggplot2)
require(sandwich)
require(MASS)

#rm(list=ls())
#-----------------------------------------------------------------------
# Read in the data and generate variables
snp.data=read.csv(file="C:/mydoc/Rutgers/Aditi/latency_snplist.csv")
#remove the summary stats that were at the bottom of the file
snp.data=snp.data[1:24,]

# create a coarser breakdown of year of diagnosis of secondary case and add to data frame
yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,48,72), include.lowest=TRUE, right=T, 
           label=c("[0,2]","(2,4]","(4,6]"))
# substitute this if you want labels in units of months
#yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,48,72), include.lowest=TRUE, right=T, 
#         label=c("[0,24]","(24,48]",">48"))

# add new breakdown of year of diagnosis to the dataframe
snp.data=data.frame(snp.data,yrgp)
rm(yrgp)
names(snp.data)
summary(snp.data)

# range of possible generation times
gentime=seq(1,340) # from 18 to 318 hours (+/- so it reaches the edges of plot)
# genome details and coverage
gensize=4411532 # assumed length of genome
coverage=0.973 #97.3% coverage of the genome
# factor to get the rate per bp*hour (assume 30.4 days/month*24 hours/day) (note I did not take the log of this quantity)
fctr =coverage*gensize*30.4*24/gentime
# print out the fctr when the generation time is 18 hours (for checking)
fctr[gentime==18]
#[1] 173986116
#--------------------------------------------------------------------------------------
# check to see if Poisson model with model based standard errors (SE) might make sense or if we need to use robust SEs
mean(snp.data$X..SNPs.different)
# 2.25
var(snp.data$X..SNPs.different)
# 10.02174 -> since variance=10.02 is much higher than the mean=2.25 we will use robust SEs

# check the Poisson model if we remove the pair with 13 SNP difference
mean(snp.data$X..SNPs.different[snp.data$X..SNPs.different<13])
# 1.782609
var(snp.data$X..SNPs.different[snp.data$X..SNPs.different<13])
# 4.996047
# the variance is still higher than the mean
# The data was checked and confirmed to be correct

#--------------------------------------------------------------------------------------
#fit the Poisson model with Months btwn diagnosis as a linear term (note: no offset)
M1=glm(formula=X..SNPs.different~Months.between.diagnosis , family="poisson", data=snp.data)

#reduced Poisson model (intercept only)
#M0=update(M1, .~. - Months.between.diagnosis)
M0=glm(formula=X..SNPs.different~1 , family="poisson", data=snp.data)

#test whether Month.between.diagnosis is significant
anova(M1,M0)
pchisq(2*(logLik(M1) - logLik(M0)), df = 1, lower.tail = FALSE)
#p=0.88

# this is the poisson estimate and 95% CI overall obtained from the reduced (intercept only) model
est.M0=c(Estimate=coef(M0),confint(M0))
# calculat robust SE and estimate with 95% CIs using them 
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

#note: the above shows p=0.90 therefore we reject the hypothesis that SNPs differ by Months.between.diagnosis

#----------------------------------------------------------------
# Estimate a rate per bp*generation
# if values of Months between diagnosis = 0 to 1 since log(0) is undefined

# need to make a new dataframe for this analysis so we are not changing the data in the original dataframe
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
#--------------------------------------------------------------------------------------
# set parameters that will be used later
# range of possible generation time
gentime=seq(1,340) # from 18 to 318 hours (+/- so it reaches the edges of plot)
gensize=4411532 # length of TB genome
coverage=0.973 #97.3% coverage of the genome
# factor to get the rate per bp*hour (assume 30.4 days/month*24 hours/day) (note I did not take the -log of this so to get rate I just divide)
fctr=coverage*gensize*30.4*24/gentime

# this is the yaxis limits for the graph
ylim=c(1e-11,1e-6)

#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis was [0,2] years/[0,24] months
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
# note this is not used further 

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

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M0to2yr=poissrate0to2yr[1]/fctr
rate.M0to2yr.LL=poissrate0to2yr[2]/fctr
rate.M0to2yr.UL=poissrate0to2yr[3]/fctr

# plot results
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\Rate0to2yr.pdf",height=6,width=4)
par(mfrow=c(1,1), mar=c(4,5,2,0)+0.1)
plot(gentime, rate.M0to2yr, 
     xlim=c(18,318), ylim=ylim, 
     xlab="Generation time (hours)", 
     ylab="Mutation Rate (mutations/bp/generation)",col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
# title 
mtext("Latency: [0, 2] Years",side=3, line= 1)

# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M0to2yr.LL,rev(rate.M0to2yr.UL)),
            col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)
# put y-axis tick marks w/o label
axis(side=2,at=c(seq(1:9)*1e-11,
      seq(1:9)*1e-10,seq(1:9)*1e-9,seq(1:9)*1e-8,seq(1:9)*1e-7,seq(1:9)*1e-6),
      labels=F,cex=.5,tcl=-.2)
# add y-axis labels only at 10^super
super=c(-11,-10,-9,-8,-7,-6)
yat=10^super
ylab=parse(text=paste(rep("10^",length(super)),super,sep=""))
axis(side=2,at=yat, label=ylab, las=2, tcl=-.3, lwd=1, adj=.5) #tcl is used to shorten tick marks
# put x-axis tick marks and labels 
xat=18+60*(0:5)
axis(side=1, at=xat, label=xat, las=1, lwd=1, adj=.5)

#add horizontal line at rate for 0-2 year interval rate for generation time of 18 hrs
abline(h=rate.M0to2yr[gentime==18],lty=2,col="grey")
dev.off()

#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis of (2,4] years/(24,48] months

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
# note this is not used further 

# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was >2-4 years
poissrate2to4yr<- c(Estimate= exp(coef(M2to4yr)), 
                    LL = exp(coef(M2to4yr) - 1.96 * std.err.M2to4yr),
                    UL = exp(coef(M2to4yr) + 1.96 * std.err.M2to4yr))
poissrate2to4yr
#Estimate.(Intercept)                   LL                   UL 
#          0.09763314           0.04120543           0.23133428 

# check to see if my estimate from glm is correct by calc SNP/months
sum(snp.data1$X..SNPs.different[snp.data$yrgp=="(2,4]"])
#33
sum(snp.data1$Months.between.diagnosis[snp.data$yrgp=="(2,4]"])
#338 months
##estimate
#33/338
#[1] 0.09763314
# matches above estimate from glm

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M2to4yr=poissrate2to4yr[1]/fctr
rate.M2to4yr.LL=poissrate2to4yr[2]/fctr
rate.M2to4yr.UL=poissrate2to4yr[3]/fctr

# plot results
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\Rate2to4yr.pdf",height=6,width=4)
par(mfrow=c(1,1), mar=c(4,5,2,0)+0.1)
plot(gentime, rate.M2to4yr, 
     xlim=c(18,318), ylim=ylim, 
     xlab="Generation time (hours)", 
     ylab="Mutation Rate (mutations/bp/generation)", col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
# title 
mtext("Latency: (2, 4] Years",side=3, line= 1)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M2to4yr.LL,rev(rate.M2to4yr.UL)),
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

dev.off()
#here
#--------------------------------------------------------------------------------------
#fit the Poisson model for HHCs whose time of diagnosis of (4+] years/(48+] months

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
# note this is not used further 

# poisson rate and robust SE based 95% CIs among those HHCs with secondary diagnosis was >4 years
poissrate4plusyr<- c(Estimate= exp(coef(M4plusyr)), 
                    LL = exp(coef(M4plusyr) - 1.96 * std.err.M4plusyr),
                    UL = exp(coef(M4plusyr) + 1.96 * std.err.M4plusyr))
poissrate4plusyr
#Estimate.(Intercept)                   LL                   UL 
#         0.021276596          0.009296854          0.048693197 
# check to see if my estimate from glm is correct by calc SNP/months
sum(snp.data1$X..SNPs.different[snp.data$yrgp=="(4,6]"])
#5
sum(snp.data1$Months.between.diagnosis[snp.data$yrgp=="(4,6]"])
#235 months
##estimate
#5/235
#[1] 0.0212766
# matches above estimate from glm

# now use the range of generation times to obtain estimates rates and robust SE based 95% CIs
rate.M4plusyr=poissrate4plusyr[1]/fctr
rate.M4plusyr.LL=poissrate4plusyr[2]/fctr
rate.M4plusyr.UL=poissrate4plusyr[3]/fctr

# plot results
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\Rate4plusyr.pdf",height=6,width=4)
par(mfrow=c(1,1), mar=c(4,5,2,0)+0.1)
plot(gentime, rate.M4plusyr, 
     xlim=c(18,318), ylim=ylim, 
     xlab="Generation time (hours)", 
     ylab="Mutation Rate (mutations/bp/generation)", col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
# title 
mtext("Latency: (4, 6] Years",side=3, line= 1)
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
legend(138,5e-11,legend=expression(paste("rate=8.36×10"^"-10")),lty=2,bty="n",cex=.75)
dev.off()

#=============================================================
# create a PDF file with all three plots side by side
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\RateAll.pdf",height=6,width=10)

par(mfrow=c(1,3), mar=c(4,5,2,0)+0.1)
#----------------------------------------------------------
# 0-2 yrs
plot(gentime, rate.M0to2yr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: [0, 2] Years",
     xlab="Generation time (hours)", 
     ylab="Mutation Rate (mutations/bp/generation)", col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
# confidence interval shaded in
polygon(x=c(gentime,rev(gentime)),y=c(rate.M0to2yr.LL,rev(rate.M0to2yr.UL)),
        col=adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)
# put y-axis tick marks w/o label
axis(side=2,at=c(seq(1:9)*1e-11,
                 seq(1:9)*1e-10,seq(1:9)*1e-9,seq(1:9)*1e-8,seq(1:9)*1e-7,seq(1:9)*1e-6),
     labels=F,cex=.5,tcl=-.2)
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

#----------------------------------------------------------
# 2-4 yrs
plot(gentime, rate.M2to4yr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: (2, 4] Years",
     xlab="Generation time (hours)", 
    # ylab="Mutation Rate (mutations/bp/generation)", 
     ylab="" ,
     col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
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

#----------------------------------------------------------
plot(gentime, rate.M4plusyr, 
     xlim=c(18,318), ylim=ylim, 
     main="Latency: (4, 6] Years",
     xlab="Generation time (hours)", 
    # ylab="Mutation Rate (mutations/bp/generation)", 
    ylab="" ,
    col="dodgerblue",
     xaxt="n", yaxt="n" ,log='y',type="l",lty=1)
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
legend(198,1e-11,legend=expression(paste("rate=8.36×10"^"-10")),lty=2,bty="n",cex=.75)
#----------------------------------------------------------
dev.off()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Create a plot using time an 18 hour generation time and generation*bp as the denominator
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
#Estimate  Robust SE     Pr(>|z|)           LL           UL
#(Intercept)              -19.64359757 0.44912212 0.000000e+00 -20.52387693 -18.76331822
#Months.between.diagnosis  -0.05320629 0.01247237 1.990696e-05  -0.07765213  -0.02876044

summary(snp.data$Months.between.diagnosis[(snp.data$Years.between.diagnosis=="0"|
            snp.data$Years.between.diagnosis=="<1")])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="1-2"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="2-3"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="3-4"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="4-5"])
summary(snp.data$Months.between.diagnosis[snp.data1$Years.between.diagnosis=="5-6"])

# calculate empirical (non-model) based estimates of rates
# numerator: total number of SNPs in each category
snpperyr=c(sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="0"|
                                             snp.data1$Years.between.diagnosis=="<1"]),
           sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="1-2"]),
           sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="2-3"]),
           sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="3-4"]),
           sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="4-5"]),
           sum(snp.data$X..SNPs.different[snp.data$Years.between.diagnosis=="5-6"]))
# denominator: total Months at risk
tot.moperyear=c(sum(snp.data1$Months.between.diagnosis[snp.data$Years.between.diagnosis=="0"|
                                                         snp.data$Years.between.diagnosis=="<1"]),
                sum(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="1-2"]),
                sum(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="2-3"]),
                sum(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="3-4"]),
                sum(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="4-5"]),
                sum(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="5-6"]))

# mean months at risk for x-axis 
mean.moperyear=c(mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="0"|
                                                          snp.data$Years.between.diagnosis=="<1"]),
                 mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="1-2"]),
                 mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="2-3"]),
                 mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="3-4"]),
                 mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="4-5"]),
                 mean(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="5-6"]))

# create plot with empirical and model based estimates
pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\RatevsTime.pdf", height=6, width=5)

# points for the empirical estimates with y-axis on log10 scale
plot(mean.moperyear,log10(snpperyr/tot.moperyear/fctr[gentime==18]), pch=19,
     xlab="Months between Index Case and HHC TB Diagnosis", 
     ylab="Mutation Rate (mutations/bp/generation)",
     yaxt="n" , ylim=c(-10.1, -8.5))
# add y-axis labels
super=c(-10,-9,-8)
yval=10^super
yat=log10(yval)
# put y-axis tick marks w/o label
axis(side=2,at=log10(c(seq(1:9)*1e-10,seq(1:9)*1e-9,seq(1:9)*1e-8)),
     labels=F,tcl=-.2)
axis(side=2,at=log10(c(5e-10,5e-9,5e-8)),
     labels=F,tcl=-.5)

ylab=parse(text = paste("10^",   super,  sep=""))

axis(side=2,at=super, label=ylab, las=2,tcl=-.7, lwd=1,adj=.5) #tcl is used to shorten tick marks
# model based estimates 
lines(c(.001,66) , 
      log10(exp(coef(M10)[1]*rep(1,2)+coef(M10)[2]*c(.001,66))),lty=2, col="grey")

#lines(snp.data1$Months.between.diagnosis , 
#        log10(exp(coef(M10)[1]*rep(1,length(snp.data1$Months.between.diagnosis))+coef(M10)[2]*snp.data1$Months.between.diagnosis )))

dev.off()        

#=============================================================
#NOT PART OF PLOT-JUST CHECKING THAT NO ERRORS
# checking the rate for generation time of 18 hours
16/(110*gensize*coverage*24*30.4/18)
#8.360124e-10

rate.M0to2yr[gentime==18]
#[1] 8.360124e-10

M0to2yrPr=glm(formula=X..SNPs.different~1 , offset= log(Months.between.diagnosis*coverage*gensize*24*30.4/18),
            family="poisson", data=snp.data1, subset=yrgp=="[0,2]")
exp(coef(M0to2yrPr))
#8.360124e-10 
# note this matches rate.M0to2yr[gentime==18]!
