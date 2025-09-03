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

#circle plot
p1 + scale_fill_gradient(low="grey90",high = "red")

#bar plot
ggVennDiagram(x, force_upset = TRUE, order.set.by = "name", order.intersect.by = "none")

##export the overlapped data
venn_data <- process_data(Venn(x))
View(venn_data)
# Get the list of all intersections
intersections <- venn_data[["regionLabel"]][["item"]]
head(intersections)  # Shows which combination corresponds to each element
intersections[["4"&"5"&"6"]]
write.csv(intersections[["4&5&6"]], "overlapped_genes_venn.csv", row.names = FALSE)


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
rows <- c("HOXA5","HOXA6",
          "HOXA7","HOXB2","HOXD4",
          "SHOX","SHOX2") 

# (optional) gene set 2
rows <- c("ASIP","MC1R","MITF","SOX10","SLC45A2") 

# (optional) gene set 3
rows <- c("HOXA2","PAX1","PAX6","PITX1","SATB1") 

# (optional) gene set 4
rows <- c("MFSD12","TMEM163","TRPM1","KIT","GPR143","TYR") 

# (optional) gene set 5
rows <- c("ATRN","ATRNL1","MLANA","TRPM1","SYT4","DMRT1") 

# extracting data frame rows 
input.filt <- input [input$gene_id %in% rows, ]  

head(input.filt)

#melt
melt_input.filt <- melt(input.filt, id = c("gene_id")) 
head(melt_input.filt)

melt_input.filt$variable<-gsub("_.*","",melt_input.filt$variable)

colnames(melt_input.filt) <- c("gene_id", "variable","TPM")

#level and order the names
#(optional) for gene set 1
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("HOXA5","HOXA6",
                                                                   "HOXA7","HOXB2","HOXD4",
                                                                   "SHOX","SHOX2"))
#(optional) for gene set 2
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("ASIP","MC1R","MITF","SOX10","SLC45A2"))

#(optional) for gene set 3
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("HOXA2","PAX1","PAX6","PITX1","SATB1"))

#(optional) for gene set 4
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("MFSD12","TMEM163","TRPM1","KIT","GPR143","TYR"))

#(optional) for gene set 5
melt_input.filt$gene_id <- factor(melt_input.filt$gene_id,levels=c("ATRN","ATRNL1","MLANA","TRPM1","SYT4","DMRT1"))


melt_input.filt$variable <- factor(melt_input.filt$variable,levels=c("MCR","MCW","MCB",
                                                                     "FCG","MSG","FSG"))
#specify colors
group.colors <- c(MCR="darkorange2", MCW="yellow", MCB="grey20",FCG="grey60",MSG="grey100",FSG="grey60")

#make plot
#(optional) for gene set 1
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 7,strip.position="right")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic(base_size = 16)+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")+
  theme(axis.title.x=element_blank())


#(optional) for gene set 2
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 5,strip.position="right",scales="free_y")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic(base_size = 16)+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")+
  theme(axis.title.x=element_blank())


#(optional) for gene set 3
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 5,strip.position="left",scales="free_y")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic(base_size = 16)+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")+
  theme(axis.title.x=element_blank())

#(optional) for gene set 4
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, nrow = 2 ,strip.position="right",scales="free_y")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic(base_size = 16)+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")+
  theme(axis.title.x=element_blank())

#(optional) for gene set 5
ggplot(melt_input.filt, aes(x=variable, y=TPM, fill=variable)) + 
  geom_dotplot(method="histodot",binaxis='y', stackdir='center',dotsize=2,stackratio=0.8)+
  facet_wrap(~gene_id, ncol = 6,strip.position="left",scales="free_y")+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=1, color="red")+
  theme_classic(base_size = 16)+theme(axis.text.x = element_text(angle = 45,vjust=0.5),legend.position="none")+
  theme(axis.title.x=element_blank())
