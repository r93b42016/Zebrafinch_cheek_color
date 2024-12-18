
#DESeq2_volcano
# All the libraries needed

setwd("/home/miller/Desktop/GW_finch_cheek_analysis/DEG_analysis")
getwd()

library("DESeq2")
library("pheatmap")
library("RColorBrewer")


#input "gene_count_stringtie_-e_bTaeGut1.4.pri_ZFcolor.csv"
input_all <- read.table(file.choose(), header=TRUE, row.names="gene_id", sep=",")
head(input_all)

#remove row sum = 0
input = input_all[ rowSums(input_all)!=0, ]
head(input)

###choose 1 below_color project
##all_lib
condition <- factor(c("FCG","MCR","MSG","FCG","FSG","MCR","MSG","FSG","MCR","MCB","MCB","MCW","MCW"))

##color 
#MCR vs FCG
input<-input_exp[,c(2,6,9,1,4)]
condition <- factor(c("R","R","R","G","G"))
head(input)
#MCR vs MCB
input<-input_exp[,c(2,6,9,10,11)]
condition <- factor(c("R","R","R","B","B"))
head(input)
#MCR vs MCW
input<-input_exp[,c(2,6,9,12,13)]
condition <- factor(c("R","R","R","W","W"))
head(input)
#FCG vs MCB
input<-input_exp[,c(1,4,10,11)]
condition <- factor(c("G","G","B","B"))
head(input)
#MCR vs MSG/FSG/FCG 
input<-input_exp[,c(2,6,9,3,7,8,5,1,4)]
condition <- factor(c("R","R","R"))
head(input)


##region
#FSG vs FCG
input<-input_exp[,c(8,5,1,4)]
condition <- factor(c("S","S","C","C"))
head(input)

#S vs C
input<-input_exp[,c(3,7,8,5,2,6,9,1,4)]
condition <- factor(c("S","S","S","S","C","C","C","C","C"))
head(input)

##gender
#FSG vs MSG
input<-input_exp[,c(8,5,3,7)]
condition <- factor(c("F","F","M","M"))
head(input)


#####start#####
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(input), condition))

#dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeqDataSetFromMatrix(countData=input, colData=coldata, design=~ condition)

dds

# Run the DESeq pipeline
dds <- DESeq(dds)

##to output dds (optional)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="xxx.csv")

## Get differential expression results
res <- results(dds, contrast=c("condition","R","G"))
summary(res)

#remove any rows with NA
res <- res[complete.cases(res),]  
summary(res)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "gene_id"
head(resdata)

## Write results
write.table(resdata, file="xxx.csv",sep=",", col.names = T, row.names = F,quote = FALSE)

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot")

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

#to output rld (optional)
rld_assay <- assay(rlog(dds, blind=FALSE))
write.table(rld_assay, file="xxx.csv",sep=",")


#plotPCA_DESeq2
plotPCA(rld, intgroup = "condition",
        ntop = 1000, returnData = FALSE)



### Volcano plot with "significant" genes labeled ###

library(ggplot2)
library("RColorBrewer")
library(ggrepel)
library(dplyr)
##https://www.rpubs.com/Knight95/volcano

#load DESeq2 opt 
res <- read.table(file.choose(), header=TRUE, sep=",")
head(res)

#define expression
res$expression = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, 
                      ifelse(res$log2FoldChange> 1 ,'S_up','C_up'),
                      'None')



#make plot
res$delabel <- ifelse(res$gene_id %in% c("SHOX", "HOXA6","HOXB2", "HOXA2",
                                           "PAX1","HOXD4","PAX6","HOXA5",
                                           "SATB1","HOXA7","PITX1","SHOX2"), res$gene_id, NA)
head(res)
ggplot(data = res, 
       aes(x = log2FoldChange, 
           y = -log10(padj), 
           colour=expression,
           label = delabel)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("red","grey", "blue"))+
  xlim(c(-10, 10)) +
  ylim(c(0,25)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="grey",lwd=0.8) +
  labs(x="log2FoldChange",
       y="-log10 (padj)",
       title="cheeks vs scalps")  +
  theme_bw()+
  theme(axis.text = element_text(
    color="black", 
    size=20))+
  geom_text_repel(max.overlaps = Inf, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
