### 03-18-2018
#script for TG macbook
#R version 3.3.1

#EdgeR analysis of Albert's alignment and rRNA depletion for PA14 biogeo RNAseq

require(dplyr)
require(xlsx)

counts <- read.csv("~/GoogleDrive/Biogeography/data/processed/counts_annotated.csv")

setwd("/Users/Tara/GoogleDrive/Stenotrophomonas/data/processed/edgeR")
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages("reshape2")
require("edgeR")
require(reshape2) 

counts.matrix <- counts
row.names(counts.matrix) <- counts$Identifier
counts.matrix <- counts.matrix[,5:21]
regions <- c("asm.middle", "asm.outer", 
            "asm.inner", "asm.middle", "asm.outer", 
            "asm.inner", "asm.middle", "asm.outer", 
            "min.inner", "min.middle", "min.outer",
            "min.inner", "min.middle", "min.outer",
            "min.inner", "min.middle", "min.outer")
regions2 <- c("middle", "outer", 
             "inner", "middle", "outer", 
             "inner", "middle", "outer", 
             "inner", "middle", "outer", 
             "inner", "middle", "outer", 
             "inner", "middle", "outer") 
             
media.type <- regions <- c("asm", "asm", 
                           "asm", "asm", "asm", 
                           "asm", "asm", "asm", 
                           "min", "min", "min",
                           "min", "min", "min",
                           "min", "min", "min")
Plates <- c("Plate1", "Plate1", 
            "Plate2", "Plate2", "Plate2",
            "Plate3", "Plate3", "Plate3",
            "Plate4", "Plate4", "Plate4",
            "Plate5", "Plate5", "Plate5",
            "Plate6", "Plate6", "Plate6")
samples <- cbind(regions, colnames(counts.matrix))

#edgeR stores data in a simple list-based data object called a DGEList. 
counts.matrix <- as.matrix(counts.matrix)
d <- DGEList(counts=counts.matrix,group=regions)
cpm<-(cpm(d)) # calculate counts per million 

##### PCOA
require(vegan)
# bray-curtis distanace indices with vegan vegdist
# permanova on distance matrix with vegan adonist
# pcoa on distance matrix with ape 
distance.matrix<-vegdist(t(cpm), method="bray")
sample.names <- colnames(cpm)
sample.info <- as.data.frame(cbind(media.type, regions2))
sample.info<-cbind(sample.info, Plates)
adonis(distance.matrix ~media.type+Plates+regions2, data=sample.info, permutations=999)
require(ape)
distance.matrix.coordinates<-pcoa(distance.matrix, correction="none", rn=NULL)
biplot(distance.matrix.coordinates)
distance.matrix.axes<-distance.matrix.coordinates$vector
distance.matrix.axes<-data.frame(distance.matrix.axes[,1:3])
require(ggplot2)
plot <- ggplot()
COLOR.PANEL <- c(
  asm.inner= "#b3cde3",
  asm.middle="#8c96c6",
  asm.outer="#88419d",
  min.inner="#fed98e",
  min.middle="#fe9929",
  min.outer ="#cc4c02"
)

COLOR.PANEL2 <- c(
  inner="#17202A",
  middle="#2C3E50",
  outer="#ABB2B9"
)


distance.matrix.axes$regions <- regions
distance.matrix.axes<-data.frame(distance.matrix.axes)
distance.matrix.axes$Name <- c("Plate1", "Plate1", 
                               "Plate2", "Plate2", "Plate2",
                              "Plate3", "Plate3", "Plate3",
                              "Plate4", "Plate4", "Plate4",
                              "Plate5", "Plate5", "Plate5",
                              "Plate6", "Plate6", "Plate6")
plot <- ggplot()
plot + geom_point(data=distance.matrix.axes,
                  aes(x=Axis.1, y=Axis.2, color=regions), size=4)+
  scale_color_manual(values = COLOR.PANEL)+
  theme(
    plot.background=element_rect(fill="black"),
    panel.background = element_rect(fill="black"),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.line=element_line(color = "white"),
    axis.text = element_text(color="white", size=14),
    axis.title=element_text(color="white", size=18),
    legend.background = element_rect(fill=
                                      "black"),
    legend.key = element_rect(fill="black"),
    legend.text = element_text(color="white")
  )+
  xlab("Axis 1 (22%)")+
  ylab("Axis 2 (16%)")+
  xlim(-0.3,0.4)+
geom_text(data=distance.matrix.axes,
                  aes(x=Axis.1+0.02,y=Axis.2+0.02, label=Name), color="White", size=5)

total.counts <- colSums(counts.matrix)
library.size <- cbind(total.counts, colnames(counts.matrix))
library.size<-data.frame(library.size)
library.size$Plate <- Plates
plot + geom_point(data=library.size, aes(
  x=Plate, 
                                y=as.numeric(as.character(total.counts)), color=regions), size=4)+
  scale_color_manual(values=COLOR.PANEL)+
  theme(
    plot.background=element_rect(fill="black"),
    panel.background = element_rect(fill="black"),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.line=element_line(color = "white"),
    axis.text = element_text(color="white", size=14),
    axis.title.y=element_text(color="white", size=18),
    axis.title.x=element_blank(),
    legend.background = element_rect(fill=
                                       "black"),
    legend.key = element_rect(fill="black"),
    legend.text = element_text(color="white")
  )+
  xlab("Plate Number")+
  ylab("Number of aligned reads")

library.size$total.counts <- sort(library.size$total.counts, decreasing=TRUE)
plot + geom_col(data=library.size, aes(x=V2, y=total.counts))


                 


total_gene<-apply(d$counts, 2, sum) 
keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep,]
dim(d)
#Recalculate library sizes of the DGEList object
d$samples$lib.size <- colSums(d$counts)
d$samples
#calcNormFactors" function normalizes for RNA composition
#it finds a set of scaling factors for the library sizes (i.e. number of reads in a sample)
#that minimize the log-fold changes between the samples for most genes 
d.norm<- calcNormFactors(d)
#plot MDS - multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples. 
plotMDS<-plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottom", as.character(unique(d$samples$group)), col=1:3, pch=20)
#create a design matrix using model.matrix function
samples <- data.frame(samples)
fac <- paste(samples$regions, sep=".")
fac <- factor(fac)
design <- model.matrix (~0+fac)
colnames(design) <- levels(fac)

#estimate dispersion
d.norm<- estimateGLMCommonDisp(d.norm, design)
d.norm <- estimateGLMTrendedDisp(d.norm, design)
d.norm <- estimateGLMTagwiseDisp(d.norm, design)

#mean variance plot
#grey = raw variance
#light blue = tagwise dispersion varaince
#dark blue = common dispersion
#black line = poission variance

meanVarPlot <- plotMeanVar( d.norm, show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )


#BCV plot
#the BCV plot shows the common trended and genewise dispersions
#as a function of average logCPM
#plotBCV(y)

#fit genewise negative binomial general linear model
fit <- glmFit(d.norm, design)
#the likelihood ratio stats are computed for the comparison of interest
# design matrix:
  # asm.inner, asm.middle, asm.outer, min.inner, min.middle, min.outer

'''
asm.inner.vs.min.inner <- glmLRT(fit, contrast=c(1,0,0,-1,0,0))
asm.inner.vs.min.inner.de <- decideTestsDGE(asm.inner.vs.min.inner , p=0.05, adjust="BH")
de.asm.inner.min.inner.tags<-rownames(d)[as.logical(asm.inner.vs.min.inner.de)]
plotSmear(asm.inner.vs.min.inner, de.tags = de.asm.inner.min.inner.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in asm R1 compared to minA R1", cex.main=.8)

'''



#Region Comparisons for minA
min.inner.vs.min.middle <- glmLRT(fit, contrast=c(0,0,0,-1,1,0))


min.inner.vs.min.middle.de <- decideTestsDGE(min.inner.vs.min.middle , p=0.05, adjust="BH")
min.inner.vs.min.middle.tags<-rownames(d)[as.logical(min.inner.vs.min.middle.de)]
plotSmear(min.inner.vs.min.middle, de.tags = min.inner.vs.min.middle.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in min R1 compared to min R2", cex.main=.8)
setwd("~/GoogleDrive/Biogeography/output/tables/edgeR_pairwise")
min.inner_min.middle<- topTags(min.inner.vs.min.middle, adjust="BH", n=Inf)
write.table(x=min.inner_min.middle, file = "min.inner_min.middle.txt", sep ="\t", quote = FALSE)

min.inner.vs.min.outer <- glmLRT(fit, contrast=c(0,0,0,-1,0,1))
min.inner.vs.min.outer.de <- decideTestsDGE(min.inner.vs.min.outer , p=0.05, adjust="BH")
min.inner.vs.min.outer.tags<-rownames(d)[as.logical(min.inner.vs.min.outer.de)]
plotSmear(min.inner.vs.min.outer, de.tags = min.inner.vs.min.outer.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in min R1 compared to min R2", cex.main=.8)
setwd("~/GoogleDrive/Biogeography/output/tables/edgeR_pairwise")
min.inner_min.outer<- topTags(min.inner.vs.min.outer, adjust="BH", n=Inf)
write.table(x=min.inner_min.outer, file = "min.inner_min.outer.txt", sep ="\t", quote = FALSE)


#Region Comparisons for asm
asm.inner.vs.asm.middle <- glmLRT(fit, contrast=c(-1,1,0,0,0,0))
asm.inner.vs.asm.middle.de <- decideTestsDGE(asm.inner.vs.asm.middle , p=0.05, adjust="BH")
asm.inner.vs.asm.middle.tags<-rownames(d)[as.logical(asm.inner.vs.asm.middle.de)]
plotSmear(asm.inner.vs.asm.middle, de.tags = asm.inner.vs.asm.middle.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in asm R1 compared to asm R2", cex.main=.8)
setwd("~/GoogleDrive/Biogeography/output/tables/edgeR_pairwise")
asm.inner_asm.middle<- topTags(asm.inner.vs.asm.middle, adjust="BH", n=Inf)
write.table(x=asm.inner_asm.middle, file = "asm.inner_asm.middle.txt", sep ="\t", quote = FALSE)

asm.inner.vs.asm.outer <- glmLRT(fit, contrast=c(-1,0,1,0,0,0))
asm.inner.vs.asm.outer.de <- decideTestsDGE(asm.inner.vs.asm.outer , p=0.05, adjust="BH")
asm.inner.vs.asm.outer.tags<-rownames(d)[as.logical(asm.inner.vs.asm.outer.de)]
plotSmear(asm.inner.vs.asm.outer, de.tags = asm.inner.vs.asm.outer.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in asm R1 compared to asm R2", cex.main=.8)
setwd("~/GoogleDrive/Biogeography/output/tables/edgeR_pairwise")
asm.inner_asm.outer<- topTags(asm.inner.vs.asm.outer, adjust="BH", n=Inf)
write.table(x=asm.inner_asm.outer, file = "asm.inner_asm.outer.txt", sep ="\t", quote = FALSE)

