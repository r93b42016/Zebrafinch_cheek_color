library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

setwd("C:\\Users\\r93b4\\Desktop\\Collaborations\\GW_collaboration\\ZF_sex_feather\\mapto_ZF_PATWV2\\figures\\fig2_color")
getwd()

#order the input according to percentage positin first!!
charts.data <- read.table(file.choose(), head=T, sep=",")
head(charts.data)

ordered_status <- factor(charts.data$status, level = c('up','non_DEGs','no_expressed' ))

#basic
p1 <- ggplot() + geom_bar(aes(y = number, x = category, fill = ordered_status), data = charts.data, stat="identity")
p1

#Adjusting color palette, axis, theme, remove legend tilte

fill <- c("#E1B378","#5F9EA0", "#808080")
p2 <- p1 + scale_fill_manual(values=fill)+theme_bw()+
  theme(legend.title = element_blank(), axis.title.x=element_blank(), 
        legend.position="bottom",
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y =element_text(colour="black", size = 12),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p2
