cl.out.path <- "\\\\ADD/YOUR/PATH"
setwd(cl.out.path)

library(scrattch.hicat)
library(Matrix)
library(scrattch.io)
library(mfishtools)
library(matrixStats)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scrattch.vis)
library(RANN)
library(WGCNA)
library(gridExtra)
library(data.table)
library(R.utils)
library(Seurat)
library(cowplot)

options(stringsAsFactors=FALSE)

############################################################################
## Read in 10x data and additional information

load("data/hippocampus_10x.RData")
norm.dat.10X    <- logCPM(datExpr)
row.names(anno) <- colnames(norm.dat.10X)
samp.dat.10X    <- anno

exclude.genes   <- scan("data/exclusiongenes_mito_sex_tissue.txt", what="character")
use.genes       <- sort(setdiff(row.names(norm.dat.10X), exclude.genes))
save(use.genes, exclude.genes, file="use.genes.rda")

ionChannels     <- scan("data/ion_channels.txt",what="character")
epState         <- setNames(c("WG1","WG1","WG4","WG4"), c("H16.06.008","H17.06.015","H16.06.009","H16.06.010")) 


############################################################################
## Assign classes using k-means clustering on all cells based on marker genes from HBA and Cell types database

# Read in predetermined genes from Allen Brain Map
highDG   <- intersect(scan("data/topDG_HBA.txt",what="character",sep="\n"),rownames(norm.dat.10X))[1:100]
highGlut <- intersect(scan("data/topGlut_CT.txt",what="character",sep="\n"),rownames(norm.dat.10X))
highGABA <- intersect(scan("data/topGABA_CT.txt",what="character",sep="\n"),rownames(norm.dat.10X))
highGlia <- intersect(scan("data/topGlia_CT.txt",what="character",sep="\n"),rownames(norm.dat.10X))
dexGenes <- unique(c(highDG,highGlut,highGABA,highGlia)) 

# Run K-means clustering on cells with enough genes detected
set.seed(1)  # For reproducibility
kpCells  <- samp.dat.10X$gene.counts.0 > 1000 
dat      <- t(norm.dat.10X[dexGenes,kpCells])
rcl      <- kmeans(as.matrix(as.data.frame(dat)),4)
rcl      <- kmeans(as.matrix(as.data.frame(dat)),4)
rcl      <- kmeans(as.matrix(as.data.frame(dat)),4)
# NOTE: this above line of code needs to be run at least 3 times to reproduce the results
# There is no reason why this should be the case, but it is.
cluster  <- norm.dat.10X[1,]*0
cluster[kpCells] <- rcl$cluster

# Assign clusters to appropriate classes based on marker gene expression
outVal   <- NULL
for(i in 1:4) {
  outVal <- cbind(outVal,Matrix::rowMeans(norm.dat.10X[c("PROX1","GAD1","SLC17A7","SLC1A3"),cluster==i]>0))
}
val <- 1
for (i in 1:4) {
  v   <- which.max(outVal[i,])
  val <- c(val, v+1)
  outVal[,v] = 0
}
classes  <- c("Low quality", "Dentate Gyrus", "Inhibitory", "Excitatory", "Non-neuronal")[order(val)]
classes  <- classes[cluster+1]    
samp.dat.10X$class_label <- classes
print(table(classes[kpCells]))
#Dentate Gyrus    Excitatory    Inhibitory  Non-neuronal 
#        12011          6386          5302          3543 # These are the values you should get.


############################################################################
## Plot relevant ion channels - This is figure 4B
ionch   <- c("KCNJ2","KCNA1","KCNB1","KCNC4","KCND2","KCNQ2","CACNA1D","CACNA1B","CACNA1H","KCNMA1","KCNN2","SCN1A","HCN1")
ionDat  <- 2^as.matrix(t(norm.dat.10X[ionch,]))-1  # Convert to linear CPM
plotDat <- data.frame(sample_name = colnames(norm.dat.10X),ionDat)  
annoDat <- data.frame(sample_name = colnames(norm.dat.10X),cluster = samp.dat.10X$class)
annoDat <- annotate_cat(annoDat,"cluster")

p <- sample_fire_plot(plotDat,annoDat,ionch,"cluster",log_scale=TRUE, max_width=25,label_height=15)
p
ggsave("class_markers_new.pdf",p,device="pdf",height=5,width=4)


############################################################################
## Subset only the cells in dentate gyrus
kpDG   <- samp.dat.10X$class_label=="Dentate Gyrus"
datDG  <- as.matrix(norm.dat.10X[,kpDG])
annoDG <- samp.dat.10X[kpDG,]

# Reorder donor by WG
annoDG$Donor = factor(annoDG$Donor, levels=c("H16.06.008","H17.06.015","H16.06.009","H16.06.010"))
annoDG <- annotate_factor(annoDG,"Donor")

print(table(annoDG$Donor_label))
#  H16.06.008 H16.06.009 H16.06.010 H17.06.015 
#        5035       1398       1853       3725

print(table(annoDG$Donor_label)/table(samp.dat.10X$Donor)[names(table(annoDG$Donor_label))])
#  H16.06.008 H16.06.009 H16.06.010 H17.06.015 
#   0.5678998  0.2992935  0.2318569  0.4848367
# There is a significantly higher fraction of DG cells in the WG1 cases vs. the WG4 cases.

kpGn <- rowSums(datDG>1) >= (dim(datDG)[2]/20)  # Keep genes expressed in at least 5% of DG cells
useGenes <- sort(rownames(datDG)[kpGn])
datDG <- datDG[useGenes,]

fwrite(signif(datDG,3), file = "DG_data_new.csv", row.names=TRUE)
file.remove('DG_data_new.csv.gz')
gzip('DG_data_new.csv',destname='DG_data_new.csv.gz')


############################################################################
## Compute statistics per donor and find potentially significant genes and output
cl <- setNames(annoDG$Donor_label,colnames(datDG))
meansDG <- get_cl_means(datDG,cl)
e1 <- epState[colnames(meansDG)] == "WG1"
e4 <- epState[colnames(meansDG)] == "WG4"
dexDG <- rowMeans(meansDG[,e4]) - rowMeans(meansDG[,e1])
sigDG <- dexDG / (1+apply(meansDG[,e4],1,sd)+apply(meansDG[,e1],1,sd))

stats <- signif(data.frame(meansDG,dexDG,sigDG)[order(sigDG),],3)
colnames(stats) <- c(paste(epState[colnames(meansDG)],names(epState[colnames(meansDG)])),"Fold Change", "Score")

fwrite(stats, file = "DG_stats_new.csv", row.names=TRUE)


############################################################################
## Plot all the ion channels

# The first three genes in this plot is part of Figure 4E
ionChannels <- c("KCNMA1","CACNA1B","KCNJ2",intersect(names(sort(sigDG)), ionChannels))

plotDat <- data.frame(donor = as.factor(paste(epState[annoDG$Donor_label],annoDG$Donor_label)),
                      WG = epState[annoDG$Donor_label])

pdf("ionChannels_new.pdf", onefile = TRUE)
for (i in ionChannels){
  plotDat$log2CPM_UMI <- datDG[i,]
  main <- paste0(i,"; mean(WG4-WG1) = ",signif(dexDG[i],3),"; score = ",signif(sigDG[i],3))
  p <- ggplot(plotDat, aes(x=donor, y=log2CPM_UMI, fill=WG)) + 
	   geom_boxplot(width=0.1, fill="white") +
	   labs(title=main) + 
	   ylim(0,15)
  grid.arrange(p)
}
dev.off()	 


############################################################################
## Plot the best ion channels

# This plot is part of Figure 1B
ionch2 <- c(intersect(rownames(stats),ionChannels)[1:10],
            intersect(rownames(stats)[dim(stats)[1]:1],ionChannels)[10:1])
ionDG  <- 2^as.matrix(t(datDG[ionch2,]))-1  # Convert to linear CPM
plotDG <- data.frame(sample_name = colnames(datDG),ionDG)  
annoDG$sample_name = colnames(datDG)
annoDG <- annoDG[,c("sample_name",setdiff(colnames(annoDG),"sample_name"))]

p = sample_fire_plot(plotDG,annoDG,ionch2,"Donor",log_scale=TRUE, max_width=35,label_height=15)
p
ggsave("ionChannel_DG_markers_new.pdf",p,device="pdf",height=5,width=2.5)


############################################################################
## Plot the best neuron projection genes

# This plot is part of Figure 1B
npGenes <- scan("data/neuron_projection_genes.txt", what="character")
np      <- c(intersect(rownames(stats),npGenes)[1:10],
            intersect(rownames(stats)[dim(stats)[1]:1],npGenes)[10:1])
npDG    <- 2^as.matrix(t(datDG[np,]))-1  # Convert to linear CPM
plotDG  <- data.frame(sample_name = colnames(datDG),npDG)  

p = sample_fire_plot(plotDG,annoDG,np,"Donor",log_scale=TRUE, max_width=35,label_height=15)
p
ggsave("neuronProjection_DG_markers_new.pdf",p,device="pdf",height=5,width=2.5)


############################################################################
## Compare WG4/WG1 DEX genes to other combos.

dexR1 <- rowMeans(meansDG[,c(3,4)]) - rowMeans(meansDG[,c(1,2)])
dexR2 <- rowMeans(meansDG[,c(2,4)]) - rowMeans(meansDG[,c(1,3)])

out_dex <- data.frame(
    real_dex  = c(sum(dexDG >= 1),sum(dexDG <= -1)),
	rand_dex1 = c(sum(dexR1 >= 1),sum(dexR1 <= -1)),
	rand_dex2 = c(sum(dexR2 >= 1),sum(dexR2 <= -1)))
rownames(out_dex) <- c("Higher in WG4", "Lower in WG4")
print(out_dex)
#               real_dex rand_dex1 rand_dex2
# Higher in WG4      397       634       303
# Lower in WG4       341       176       289
# The number of genes is no different from what we'd expect by chance.  What if we put all the cells on a PC vector?


############################################################################
## Run Seurat on the data and output plots for PCs and UMAP

human_metadata <- data.frame(donor = annoDG$Donor_label,
                             state = epState[annoDG$Donor_label])
rownames(human_metadata) <- colnames(datDG)
human <- CreateSeuratObject(datDG, meta.data = human_metadata)
human <- NormalizeData(human, verbose = FALSE)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000, verbose = FALSE) # 750
human <- ScaleData(object = human, verbose = FALSE)
human <- RunPCA(object = human, npcs = 30, verbose = FALSE)
human <- RunUMAP(object = human, reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(object = human, group.by = "donor", reduction = "pca", do.return = TRUE, 
      pt.size = 0.5, label=TRUE) + NoLegend() 
p2 <- DimPlot(object = human, group.by = "state", reduction = "pca", do.return = TRUE, 
      pt.size = 0.5, label=TRUE, cols=c("blue","red")) + NoLegend()
p3 <- DimPlot(object = human, group.by = "donor", reduction = "umap", do.return = TRUE, 
      pt.size = 0.5, label=TRUE) + NoLegend() 
p4 <- DimPlot(object = human, group.by = "state", reduction = "umap", do.return = TRUE, 
      pt.size = 0.5, label=TRUE, cols=c("blue","red")) + NoLegend()
plot_grid(p1, p2, p3, p4, ncol = 2)
ggsave("umap_pca_plots.pdf",device="pdf",height=9,width=9)

# This plot is part of Figure 1B
DimPlot(object = human, group.by = "state", reduction = "pca", do.return = TRUE, 
      pt.size = 0.4, label=FALSE, cols=c("blue","red")) + NoLegend()
ggsave("umap_pca_plots_PC_byState.pdf",device="pdf",height=5,width=5)

