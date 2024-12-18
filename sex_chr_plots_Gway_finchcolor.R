library(tidyverse)
library(dplyr)

setwd("/home/miller/Desktop/GW_finch_cheek_analysis/Fig. 5_sex_chr_0724")
getwd()

#load S4 file MCRvsFCG
input_1 <- read.table(file.choose(), header=T, sep=",")
head(input_1)

#load "bTaeGut1.4.pri_genome_stat.csv"
chr_gene_stat <- read.table(file.choose(), header=T, sep=",")
head(chr_gene_stat)

#count DEG numbers in each chr
input_1_sum <- input_1 %>% 
  group_by(chr) %>% 
  summarise(chr_sum = n())
            
head(input_1_sum)            
            
#divide DEG num per chr by toal gene num per chr
merge<-merge(input_1_sum,chr_gene_stat,by="chr")
head(merge)
            
cheek<-select(merge, chr,chr_sum, Gene)%>%
  mutate(ratio=chr_sum/Gene,region="C")

head(cheek)

#load S4 file MCRvsFCG
input_2 <- read.table(file.choose(), header=T, sep=",")
head(input_2)

#count DEG numbers in each chr
input_2_sum <- input_2 %>% 
  group_by(chr) %>% 
  summarise(chr_sum = n())

head(input_2_sum)            

#divide DEG num per chr by toal gene num per chr
merge<-merge(input_2_sum,chr_gene_stat,by="chr")
head(merge)

head<-select(merge, chr,chr_sum, Gene)%>%
  mutate(ratio=chr_sum/Gene,region="S")

head(head)


#make ratio table
input <- rbind(cheek, head)
head(input)

###make boxplot with permutation test###
#define expression
input$sex_chr = ifelse(input$chr=="W"|input$chr=="Z",'Sex_chr','Autosome')
unique(input$sex_chr)
View(input)

ggplot(input,aes(sex_chr,ratio,color=region))+ 
  geom_boxplot()+
  geom_jitter()+
  scale_color_manual(values = c("S" = "grey40", "C" = "darkorange2"))+
  facet_wrap(~region, ncol = 2,scales="free_y")+
  theme_classic()+theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(color="black",size=14))

#permutation test
library(lmPerm)
summary(lmp(ratio~sex_chr,data=input))


###point_plot_DEG ratio###
p<-ggplot(data=input, mapping=aes(x=chr, y=ratio, color=region))+
  geom_point(size=3)+
  geom_line(aes(group = region)) + 
  scale_color_manual(values = c("H" = "grey40", "C" = "darkorange2"))+
  labs(title="DEG ratio",x="Chromosme", y = "DEGs/Total Genes Per Chr.")

p+facet_grid(region ~ .)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5))

###vln_plot###
# geom_boxplot proposes several arguments to custom appearance
#first
install.packages("ggbreak")
library(ggrepel)
library(ggbreak)

head(input_1)

cheek<-select(input_1, gene_id,chr,log2FoldChange)%>%
  mutate(region="C")

head(cheek)

head<-select(input_2, gene_id,chr,log2FoldChange)%>%
  mutate(region="H")

head(head)


input <- rbind(cheek, head)

head(input)
unique(input$chr)

sub<-input %>% filter (chr=="Z"|chr=="W") %>% 
  filter(abs(log2FoldChange)>2)%>%
  filter(!grepl('LOC', gene_id))
head(sub)
unique(sub$chr)
unique(sub$gene_id)


#https://ggrepel.slowkow.com/articles/examples
ggplot(data=input,mapping=aes(x=region, y=log2FoldChange, color=region, label=gene_id)) + 
  geom_violin()+
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  #geom_point(size=3)+
  scale_color_manual(values = c("darkorange2","grey40"))+
  #scale_y_break(c(-1, -12))+
  labs(title="DEGs on sex chromosomes")+
  #geom_text_repel(data=sub, color="black", min.segment.length = Inf)+
  geom_hline(yintercept=2, linetype="dashed", color="blue",linewidth=0.5)+
  geom_hline(yintercept=0, linetype="dashed", color="black",linewidth=0.5)+
  geom_hline(yintercept=-2, linetype="dashed", color="blue",linewidth=0.5)+
  theme_classic()

  

  
  
               
 