# Rscript for pvalue calculation based on Possion distribution
library(stats)
args = commandArgs (T)
InputBedFile = args[1]
Outfile = args[2]

data <- read.table(InputBedFile)
# This ppois function will calculate a pvalue for each site using the expected count of gene E[X] based on poisson distrubution of non-zero positions.
# Expected count stored in col 9. will be replaced with pvalue.
data$V9 <- ppois(data$V5,lambda=data$V9,lower.tail=FALSE)#Trying mean as lambda (mean of non-zero positions).

data <- data[which(data$V9 < 0.05),]
data <- data[which(data$V7 < 0.001),]

write.table(data, file=Outfile,sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

