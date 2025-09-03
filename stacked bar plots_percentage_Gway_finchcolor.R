http://t-redactyl.io/blog/2016/01/creating-plots-in-r-using-ggplot2-part-4-stacked-bar-plots.html

library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

setwd("C:\\Users\\r93b4\\Desktop\\Collaborations\\GW_collaboration\\ZF_sex_feather\\mapto_ZF_PATWV2\\figures")
getwd()

#order the input according to percentage positin first!!
charts.data <- read.table(file.choose(), head=T, sep=",")
head(charts.data)

ordered_status <- factor(charts.data$comp, level = c('MvsF_H','MvsF_C','RvsB_C','RvsW_C','GvsB_C','GvsW_C'))
ordered_chr <- factor(charts.data$chr, level = c('Z','W','others'))

#basic
p1 <- ggplot() + geom_bar(aes(y = percentage
                            , x = ordered_status, fill = ordered_chr), data = charts.data, stat="identity")
p1
#Adding data labels
p2 <- p1 + geom_text(data=charts.data, aes(x = ordered_status, y = percentage,
                                           label = paste0(percentage,"%")), size=4)
p2
#Adjusting data labels position
charts.data <- ddply(charts.data, .(percentage), transform, pos = cumsum(percentage) - (0.5 * percentage))
head(charts.data)

p3 <- p1 + geom_text(data=charts.data, aes(x = ordered_status, y = pos, label = paste0(percentage,"%")),
                     size=4)
p3

#Adjusting color palette, axis, theme, remove legend tilte

fill <- c("#5F9EA0", "#E1B378","#808080")
p4 <- p1 + scale_fill_manual(values=fill)+theme_bw()+
  theme(legend.title = element_blank(), axis.title.x=element_blank(), 
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y =element_text(colour="black", size = 12),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p4
