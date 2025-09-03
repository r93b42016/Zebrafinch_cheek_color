# Plot heatmap

install.packages("pheatmap")

# load package
library(pheatmap)
library(grid)
library("RColorBrewer")

#setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\Gway_sex_finch\\for_heatmap")
#getwd()

setwd("/home/miller/Desktop/GW_finch_cheek_analysis")
getwd()

#input "UQTPM"
input_all <- read.table(file.choose(), header=TRUE, row.names="gene_id", sep=",")
head(input_all)
summary(input_all)

input<-input_all[,c(1:9)]
head(input)

# remove read counts sum rows< X
input_exp = input[ rowSums(input)>3, ]
head(input_exp)



#for colorDEG
input<-input_all[,c(2,3,5,7,8,9,10,11,12,13,14,15,16,17)]
head(input)

#for allcompDEG
input<-input_all[,c(2,3,5,7,8,9,10,11,12,13)]
head(input)

###for_dds
#replace zero to small value
input[input < 0.00001] <- 0.000001
head(input)
#log transform (optional)
input_log10<-log10(input)
head(input_log10)
#z-score transform
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
#z-score_matrix
input_UQTPMZ <- t(apply(input_log10, 1, cal_z_score))


#remove NA from z-score matrix
row.has.na <- apply(input_UQTPMZ, 1, function(x){any(is.na(x))})
sum(row.has.na)
input_UQTPMZ <- input_UQTPMZ[!row.has.na,]
write.table(input_UQTPMZ,"whether_NA")



#basic_pheatmap
out<-pheatmap(input_exp,
              scale="row",
              show_rownames = F,
              cluster_rows = T,
              clustering_method="ward.D2"
)


##add_annotation_or_condition
annotdf<- read.table(file.choose(), header=TRUE, row.names="gene_id", sep=",")
annotdf
#add_basic_color
mycolors = list(category = c(carotenoid = "red", melanogenesis = "grey"))
mycolors = list(category = c(all = "black", sex_DEG = "Magenta", color_DEG = "green", region_DEG="brown", color_region ="orange", sex_color ="blue", sex_region ="yellow"))
mycolors

#or by import list
x <- readLines(file.choose())
selected_labels <- strsplit(x,split=",")

head(selected_labels)

head(input)
summary(input)

#cut the tree at a pre-defined tree height, and extract the gene-to-cluster assignments at that height
plot(out$tree_row)
abline(h=800, col="red", lty=2, lwd=2)

clusters<-sort(cutree(out$tree_row, k=10))
head(clusters)
tail(clusters)

#check cluster number
#make cluster table for pheatmap cutting
my_gene_col <- data.frame(clusters)
head(my_gene_col)
tail(my_gene_col)
my_gene_col$clusters <- paste("cluster", my_gene_col$clusters, sep="_")
head(my_gene_col)
tail(my_gene_col)

#only for output
head(my_gene_col)
colnames(my_gene_col)[1] <- "gene_id"
colnames(my_gene_col)[1]

colnames(my_gene_col)[2] <- "clusters"
head(my_gene_col)
write.table(my_gene_col,"Rplot_cluster_wardD2_k10_GWZF_allcomp_nomutant_DEG_rmlib_UQTPM_1112.csv", sep=",")

###make color table
#for color less than 9 (not finish yet!!!)
unique_clusters<- list(unique(my_gene_col$clusters))
mycolors2 <- list(brewer.pal(7, "Set1")[1:7])
mycolors2
head(my_gene_col$clusters)
length(mycolors2)
names(mycolors2) <- unique(my_gene_col$clusters)
head(mycolors2)

#for color more than 12
my_gene_col <- data.frame(clusters)
head(my_gene_col)
my_gene_col$clusters <- paste("cluster", my_gene_col$clusters, sep="_")
head(my_gene_col)
tail(my_gene_col)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(my_gene_col$clusters))))
mycolors <- newCols(length(unique(my_gene_col$clusters)))
names(mycolors) <- unique(my_gene_col$clusters)
mycolors <- list(clusters = mycolors)

#coloered pheatmap
heatmap_inpu<-pheatmap(input_UQTPMZ,
                       scale="row",
                       #cluster_rows = FALSE,
                       #cluster_cols = FALSE,
                       #gaps_row=c(3,71),
                       #cellheight = 4.5,
                       #cellwidth = 20,
                       #border_color=NA,
                       #fontsize=6,
                       #fontsize_row = 5,
                       #fontsize_col = 6,
                       #main="Heatmap for expressed genes",
                       annotation_row = my_gene_col,
                       cutree_rows = 10,
                       treeheight_col=15,
                       treeheight_row=20,
                       show_rownames = F,
                       annotation_colors = mycolors,
                       clustering_method="ward.D2"
)

heatmap_inpu

###only display selected row.names
add.flag(heatmap_inpu,
         kept.labels = selected_labels,
         repel.degree = 0.5)

#function behind (execute prior to above!!!)
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

