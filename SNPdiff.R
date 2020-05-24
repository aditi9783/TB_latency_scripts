# This program plots the raw data and a loess smoother
rm(list=ls)

#install.packages("sandwich")
#install.packages("MASS")
#update.packages()
require(ggplot2)
require(sandwich)
require(MASS)

snp.data=read.csv(file="C:/mydoc/Rutgers/Aditi/latency_snplist.csv")
#remove the summary stats that were at the bottom of the file
snp.data=snp.data[1:24,]

# create a coarser breakdown (compared to Years.between.diagnosis) of year of diagnosis of secondary case
yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,48,72), include.lowest=TRUE, right=T, 
           label=c("[0,2]","(2,4]","(4,6]"))
# substitute this if you want labels in units of months
#yrgp=cut(snp.data$Months.between.diagnosis, breaks=c(0,24,48,72), include.lowest=TRUE, right=T, 
#         label=c("[0,24]","(24,48]",">48"))

# add new breakdown of year of diagnosis to the dataframe
snp.data=data.frame(snp.data,yrgp)
names(snp.data)
summary(snp.data)

#----------------------------------------------------------------
# plot with vertical lines at 1 year intervals 

# first get summaries (not used for the plot)
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="<1"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="1-2"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="2-3"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="3-4"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="4-5"])
summary(snp.data$Months.between.diagnosis[snp.data$Years.between.diagnosis=="5-6"])
#note there are none at exactly 12,36 mos but some at exactly 0,24,60 months

# these are n's for the plots with 6 categories displayed
n=summary(snp.data$Years.between.diagnosis)
n=n[-1] #remove the 0 in front
neq=paste("n=",n,sep="")
lab=c(expression("[0, 1]\n"),"(1, 2]\n","(2, 3]\n","(3, 4]\n","(4, 5]\n","(5, 6]\n")
# substitute below if you want months displayed
#lab=c(expression("[0, 12]\n"),"(12, 24]\n","(24, 36]\n","(36, 48]\n","(48, 60]\n",">60\n")
lab=paste(lab,neq,sep=" ")

#plot the months btwn diagnosis and SNP data with a scatterplot smoother

pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\SNPvsMonths.pdf")
par(mar=c(6,4,5,2)+0.1)

with(snp.data, scatter.smooth(Months.between.diagnosis, X..SNPs.different,
     family="symmetric",
     xlim=c(-1,66),ylim=c(0,15), xaxt="n", yaxt="n" , pch=20,
     xlab="Months between Index Case and HHC TB Diagnosis" ,
     ylab="Number of SNPs"))

axis(side=1,at=seq(0,66,by=6))
axis(side=2,at=seq(0,14,by=1),las=1)
yloc=15 # where along y-axis to put text 
text(3,yloc,"Years*:\nPairs:",adj=c(1,NA))
#text(3,yloc,"Mos.*:\nPairs:",adj=c(1,NA))
text(8,yloc,label=lab[1])
abline(v=12,lty=2)
text(18,yloc,label=lab[2])
abline(v=24,lty=2)
text(30,yloc,label=lab[3])
abline(v=36,lty=2)
text(42,yloc,label=lab[4])
abline(v=48,lty=2)
text(54,yloc,label=lab[5])
abline(v=60,lty=2)
text(66,yloc,label=lab[6])

mtext("*() exclusive, [] inclusive of value", side=1, line=5, adj=1)
# substitute the next line if the text should be inside the plot region
#text(68,13.5,label="*() exclusive\n[] inclusive", adj=1)

dev.off()
#----------------------------------------------------------------
# plot with vertical lines at 2 year intervals

# first get summaries (not used for the plot)
summary(snp.data$Months.between.diagnosis[snp.data$yrgp=="[0,2]"])
summary(snp.data$Months.between.diagnosis[snp.data$yrgp=="(2,4]"])
summary(snp.data$Months.between.diagnosis[snp.data$yrgp=="(4,6]"])

#these are n's for the coarser grouping
ngp=summary(snp.data$yrgp)
ngpeq=paste("n=",ngp,sep="")
ngpeq
labgp=c(expression("[0, 2]\n"),"(2, 4]\n","(4, 6]\n")
#labgp=c(expression("[0, 24]\n"),"(24, 48]\n",">48\n")
labgp=paste(labgp,ngpeq,sep=" ")

#plot the months btwn diagnosis and SNP data

pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\SNPvsMonths2.pdf")
par(mar=c(6,4,5,2)+0.1)
with(snp.data, scatter.smooth(Months.between.diagnosis, X..SNPs.different,
     family="symmetric",
     xlim=c(-1,66),ylim=c(0,15), xaxt="n", yaxt="n" , pch=20,
     xlab="Months between Index Case and HHC TB Diagnosis" ,
     ylab="Number of SNPs"))

axis(side=1,at=seq(0,66,by=6))
axis(side=2,at=seq(0,14,by=1),las=1)
yloc=15 # where along y-axis to put text 
text(6,yloc,"Years*:\nPairs:",adj=c(1,NA))
#text(6,yloc,"Mos.*:\nPairs:",adj=c(1,NA))
text(12,yloc,label=labgp[1])
abline(v=24,lty=2)
text(36,yloc,label=labgp[2])
abline(v=48,lty=2)
text(60,yloc,label=labgp[3])

mtext("*() exclusive, [] inclusive of value", side=1, line=5, adj=1)
# substitute the next line if the text should be inside the plot region
#text(68,13.5,label="*() exclusive\n[] inclusive", adj=1)

dev.off()

