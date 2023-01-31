#Chi-Square test for differential PolyA usage
options("scipen"=100, "digits"=5)
library("gplots")
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(Biobase)
library(reshape2)

args = commandArgs (T)
InputBedFile = args[1]
Outfile = args[2]
NoC1=as.numeric(args[3])
NoC2=as.numeric(args[4])
OutDIR = args[5]
NoS=NoC1+NoC2

data <- read.delim(InputBedFile,header=FALSE)
#data <- read.delim(InputBedFile,header=TRUE)
head(data)
data <- data[rowSums(data[,8:(7+NoS)]) > 10,  ]
#Filter low count PA sites!

genes <- data[,NoS+8]
gene.table <- as.data.frame(table(genes))
genelist <- noquote(droplevels(gene.table[which(gene.table$Freq >1),]))
rownames(genelist) <- genelist$genes
print("Starting to Switch Testing!")
#Switch Test
tested.Genes <- NULL  
for (i in genelist$genes){ 
  nPA=genelist[i,]$Freq
  GeneTotalC1 <- sum(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  GeneTotalC2 <- sum(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  PAsC1 <- rowSums(data[which(data[,NoS+8] == i  ),8:(7+NoC1)])
  PAsC2 <- rowSums(data[which(data[,NoS+8] == i  ),(8+NoC1):(7+NoS)])
  
  for (j in 1:nPA){
    x <- round(matrix(c(PAsC1[j]/NoC1, PAsC2[j]/NoC2, (GeneTotalC1-PAsC1[j])/NoC1, (GeneTotalC2-PAsC2[j])/NoC2), byrow = TRUE, 2, 2))
    print(i)
    print(x)
    cqtest <- chisq.test(x)
    #tested.PAsites <- as.data.frame(data[which(data[,NoS+8] == i  ),c(1,2,3,4,5,6,(8+NoS),(9+NoS),(10+NoS))][j,])
    tested.PAsites <- as.data.frame(data[which(data[,NoS+8] == i  ),c(1,2,3,4,5,6,(8+NoS),(9+NoS))][j,])
    tested.PAsites$Num_Total_PA_sites_of_Gene <- nPA
    tested.PAsites$PA_MeanCond1 <- PAsC1[j]/NoC1
    tested.PAsites$PA_MeanCond2 <- PAsC2[j]/NoC2
    tested.PAsites$Rest_MeanCond1 <- (GeneTotalC1-PAsC1[j])/NoC1
    tested.PAsites$Chi2.statistics  <- cqtest$statistic
    tested.PAsites$Chi2.pvalue <- cqtest$p.val
    tested.PAsites$Fis.pvalue <- fisher.test(x)$p.value
    
    tested.Genes <- rbind(tested.Genes, tested.PAsites)
    
  }
  
}


tested.Genes$Chi.padj <- p.adjust(tested.Genes$Chi2.pvalue, n=length(tested.Genes[,1]), method="BH")
tested.Genes$Fis.padj <- p.adjust(tested.Genes$Fis.pvalue, n=length(tested.Genes[,1]), method="BH")
write.table(tested.Genes[order(tested.Genes$Chi.padj),], file=Outfile,sep="\t")
print("Switch Testing is complete!")
print("Starting to UTR Length comparison!")
#3'UTR length
data <- read.delim(paste(OutDIR,"/","3utr_CPM.bed",sep=""),header=FALSE)

rownames(data) <- data$V10
data <- data[rowSums(data[,14:(13+NoS)]) > 100,  ]
head(data,2)
genes <- data[,5]
gene.table <- as.data.frame(table(genes))
genelist <- noquote(droplevels(gene.table[which(gene.table$Freq > 1),]))
rownames(genelist) <- genelist$genes




tested.Genes <- NULL
for (i in genelist$genes){
  nPA=genelist[i,]$Freq
  print(data[which(data[,5] == i  ),5]);
  utrs <- droplevels(unique(data[which(data[,5] == i  ),5]))
  for (u in utrs){
    u.data <- data[which(data[,5] == u  ),]
    if(u.data$V6[1] == "-"){
      u.data$PAS.dist <- u.data$V3 - (u.data$V8 +(u.data$V9-u.data$V8)/2)
      u.data <- u.data[order(rowMeans(data[which(data[,5] == u  ),14:(13+NoS)]),decreasing = TRUE),]
      u.data <- u.data[1:2,]
      u.data <- u.data[order(u.data$PAS.dist),]
      u.data$Site <- c("Proximal","Distal")
    }else{
      u.data$PAS.dist <- (u.data$V8 +(u.data$V9-u.data$V8)/2)-u.data$V2
      u.data <- u.data[order(rowMeans(data[which(data[,5] == u  ),14:(13+NoS)]),decreasing = TRUE),]
      u.data <- u.data[1:2,]
      u.data <- u.data[order(u.data$PAS.dist),]
      u.data$Site <- c("Proximal","Distal")
    }
    u.data$C1mean <- rowMeans(u.data[which(u.data[,5] == u  ),14:(13+NoC1)])
    u.data$C2mean <- rowMeans(u.data[which(u.data[,5] == u  ),(14+NoC1):(13+NoS)])
    u.data$C1utrlen <- (u.data$C1mean+0.0001)*u.data$PAS.dist
    u.data$C2utrlen <- u.data$C2mean*u.data$PAS.dist
    utr.len.cor <- cor(u.data$C1utrlen,u.data$C2utrlen)
    plot.data <- rbind(round(t(u.data[,14:(13+NoC1)])), round(t(u.data[,(14+NoC1):(13+NoS)])))
    plot.data<- data.frame(c(rep("Normal",NoC1),rep("PE",NoC2)) ,plot.data)
    colnames(plot.data) <- c("Type","Proximal","Distal")
    plot.data$P2D_ratio <- plot.data$Proximal/plot.data$Distal
    plot.data <- plot.data[is.finite(plot.data$P2D_ratio),]
    u.data$P2D_Ratio_N <- median(plot.data[plot.data$Type == "Normal",]$P2D_ratio)
    u.data$P2D_Ratio_PE <- median(plot.data[plot.data$Type == "PE",]$P2D_ratio)
    tested.Genes <- rbind(tested.Genes,u.data)
  }
}
print(head(tested.Genes));
tested.Genes$P2D_N_2_P2D_PE <- log2(tested.Genes$P2D_Ratio_N/tested.Genes$P2D_Ratio_PE)
tested.Genes <- tested.Genes[order(tested.Genes$P2D_N_2_P2D_PE),]
tested.Genes <- tested.Genes[abs(tested.Genes$P2D_N_2_P2D_PE)>1.2,]
write.table(tested.Genes,file=paste(OutDIR,"/","UTRlen.txt",sep=""),sep="\t")

nms <- as.data.frame(table(droplevels(tested.Genes[,5])))

pdf(paste(OutDIR,"/","UTR_Length_Shorthening.pdf",sep=""))
for (i in nms$Var1 ){
  nPA=genelist[i,]$Freq
  utrs <- droplevels(unique(data[which(data[,5] == i  ),5]))
  for (u in utrs){
    u.data <- data[which(data[,5] == u  ),]
    if(u.data$V6[1] == "-"){
      u.data$PAS.dist <- u.data$V3 - (u.data$V8 +(u.data$V9-u.data$V8)/2)
      u.data <- u.data[order(rowMeans(data[which(data[,5] == u  ),14:(13+NoS)]),decreasing = TRUE),]
      u.data <- u.data[1:2,]
      u.data <- u.data[order(u.data$PAS.dist),]
      u.data$Site <- c("Proximal","Distal")
    }else{
      u.data$PAS.dist <- (u.data$V8 +(u.data$V9-u.data$V8)/2)-u.data$V2
      u.data <- u.data[order(rowMeans(data[which(data[,5] == u  ),14:(13+NoS)]),decreasing = TRUE),]
      u.data <- u.data[1:2,]
      u.data <- u.data[order(u.data$PAS.dist),]
      u.data$Site <- c("Proximal","Distal")
    }
    u.data$C1mean <- rowMeans(u.data[which(u.data[,5] == u  ),14:(13+NoC1)])
    u.data$C2mean <- rowMeans(u.data[which(u.data[,5] == u  ),(14+NoC1):(13+NoS)])
    u.data$C1utrlen <- (u.data$C1mean+0.0001)*u.data$PAS.dist
    u.data$C2utrlen <- u.data$C2mean*u.data$PAS.dist
    utr.len.cor <- cor(u.data$C1utrlen,u.data$C2utrlen)
    plot.data <- rbind(round(t(u.data[,14:(13+NoC1)]))+0.0001, round(t(u.data[,(14+NoC1):(13+NoS)]))+0.0001)    
    plot.data<- data.frame(c(rep("Normal",NoC1),rep("Case",NoC2)) ,plot.data)
    colnames(plot.data) <- c("Type","Proximal","Distal")
    plot.data$P2D_ratio <- plot.data$Proximal/plot.data$Distal
    plot.data <- plot.data[is.finite(plot.data$P2D_ratio),]
    u.data$P2D_Ratio_N <- median(plot.data[plot.data$Type == "Normal",]$P2D_ratio)
    u.data$P2D_Ratio_PE <- median(plot.data[plot.data$Type == "PE",]$P2D_ratio)
    print(ggplot(data=plot.data,aes(x=Proximal,y=Distal))+geom_jitter(aes(color=Type))
          +ggtitle(paste(u.data[which(u.data[,5] == u  ),4],"\n",u,"\n","log2(D2P ratio (N/C)):",round(tested.Genes[which(tested.Genes[,5] == u  ),c("P2D_N_2_P2D_PE")][1],digits = 2)  ))) 
    
  }
}
dev.off()

print("UTR length comparison is complete!")

