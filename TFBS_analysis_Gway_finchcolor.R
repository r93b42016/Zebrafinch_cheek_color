
library(ggplot2)
library("RColorBrewer")
library(ggrepel)


#enrichment_dot_plot
# plot: dot plot
setwd("/home/miller/Desktop/GW_finch_cheek_analysis/in_silico_TFBS_prediction")
getwd()

input <- read.table(file.choose(),sep=",",header=T)
head(input)
#level the gene names

input$TFs <- factor(input$TFs,levels=c("PAX1","PAX6","PITX1","SOX10"))
input$species <- factor(input$species,levels=c("V_macroura","V_chalybeata","T_guttata","L_striata","O_arfaki","P_nigricollis","M_ater","M_alba","P_himalayana","P_montanus","E_fucata","C_ornatus","S_canaria","G_gallus"))
head(input)

#define dot shape
group.shapes <- c(show=16, no=17)


p <- ggplot(input, aes(colorGene, TFs, color=TFs, size=count_q005, shape=ishere, alpha=PCC))
p + geom_point()+theme_bw()+ facet_grid(species ~ .)+
  scale_color_brewer(palette = "Set1")+
  scale_shape_manual(values=group.shapes)+
  theme(axis.text.x = element_text(color = "black",size = 14, angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.text.y = element_text(color = "black",size = 14))+ 
  theme(axis.title.y = element_blank())+
  theme(strip.text.y = element_blank())+
  guides(col = FALSE)+
  theme(legend.position="left")

