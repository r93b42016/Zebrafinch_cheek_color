###venndiagram###

library(ggVennDiagram)
library(ggplot2)

setwd("/home/miller/Desktop/GW_finch_cheek_analysis/Venn")
getwd()

x <- read.table(file.choose(), head=T, sep=",")
head(x)

p1 <- ggVennDiagram(x)
p1

# four dimension venn plot
theme_classic()
ggVennDiagram(x, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none")




### make dotplots ###

library("ggplot2")
library(dplyr)
library(reshape2) 

setwd("/home/miller/Desktop/GW_finch_cheek_analysis/UQTPM_folder")
getwd()

#load sup file1 UQTPM table
input_all <- read.table(file.choose(), header=T, sep=",",stringsAsFactors = F)
head(input_all)

input<-input_all[,c(1,4:16)]
head(input)


# getting rows  
# (optional) gene set 1
rows <- c("SHOX","HOXA6","HOXB2",
          "HOXA2","PAX1","HOXD4",
          "PAX6","HOXA5","SATB1",
          "HOXA7","PITX1","SHOX2") 

# (optional) gene set 2
rows <- c("ASIP","MC1R","MITF","SOX10","SLC45A2") 


# extracting data frame rows 
input.filt <- input [input$gene_id %in% rows, ]  

head(input.filt)

#melt
melt_input.filt <- melt(input.filt, id = c("gene_id")) 
head(melt_input.filt)
View(melt_input.filt)

melt_input.filt$variable<-gsub("_.*","",melt_input.filt$variable)
View(melt_input.filt)

colnames(melt_input.filt) <- c("gene_id", "variable","TPM")

#level and order the names
#(optional) for gene set 1
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("HOXA2","HOXA5","HOXA6",
                                                                   "HOXA7","HOXB2","HOXD4",
                                                                   "PAX1","PAX6","PITX1",
                                                                   "SATB1","SHOX","SHOX2"))
#(optional) for gene set 2
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("ASIP","MC1R","MITF","SOX10","SLC45A2"))


melt_input.filt$variable <- factor(melt_input.filt$variable,levels=c("MCR","MCW","MCB",
                                                                     "FCG","MSG","FSG"))
#specify colors
group.colors <- c(MCR="darkorange2", MCW="yellow", MCB="grey20",FCG="grey60",MSG="grey100",FSG="grey60")

#make plot
#(optional) for gene set 1
ggplot(melt_input.filt, aes(x=variable, y=log10(TPM), fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 6,strip.position="right")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_bw()+theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank())


#(optional) for gene set 2
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 1,strip.position="right",scales="free_y")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")



