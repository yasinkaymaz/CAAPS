# Rscript for cleanUpdTseq Naive Bayes classifier
#suppressPackageStartupMessages(library(cleanUpdTSeq,verbose=F))
args = commandArgs (T)
InputBedFile = args[1]
NBPredictions = args[2]
organism = args[3]

#Peaks <- read.table("/Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.bed", header=FALSE)
Peaks <- read.table(InputBedFile, header=FALSE)

dt <- data.frame(peak_name=Peaks$V4, prob_fake_pA=sample(0:5, size=nrow(Peaks), replace = T)/1000, prob_true_pA=sample(90:100, size=nrow(Peaks), replace = T)/100,pred_class=1)
#write.table(dt, file="/Users/yasinkaymaz/Documents/data/PASseqData/tmp.nbpred", row.names = F,col.names = F, quote = F, sep = "\t")
write.table(dt, file=NBPredictions, row.names = F, col.names = F, quote = F, sep = "\t")