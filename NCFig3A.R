# NCFig3A.R 
# generate Figure 3A
# raw SNP data with a loess smoothed curve

rm(list=ls())

#----------------------------------------------------------------
# read in the data and create categories
snp.data=read.csv(file="C:/mydoc/Rutgers/Aditi/latency_snplist.csv")
#remove the summary stats that were at the bottom of the file
snp.data=snp.data[1:24,]

names(snp.data)
summary(snp.data)
#----------------------------------------------------------------
# create a pdf plot

pdf(file = "C:\\mydoc\\Rutgers\\Aditi\\RAnalysis\\NCFig3A.pdf")
par(mar=c(6,6,5,2)+0.1)

with(snp.data, scatter.smooth(Months.between.diagnosis, X..SNPs.different,
     family="symmetric",
     xlim=c(-1,64),ylim=c(0,14), xaxt="n", yaxt="n" , pch=20,
     xlab="Months between Index Case and HHC TB Diagnosis" ,
     ylab="Number of SNPs"))

axis(side=1,at=seq(0,66,by=6))
axis(side=2,at=seq(0,14,by=1),las=1)

dev.off()

