# 3-19-18

# find trends in gene regulation among regions
setwd("~/GoogleDrive/Biogeography/output/tables/edgeR_pairwise/")
asm.inner_asm.middle <- read.delim("asm.inner_asm.middle.txt")
asm.inner_asm.outer <- read.delim("asm.inner_asm.outer.txt")
min.inner_min.outer <- read.delim("min.inner_min.outer.txt")
counts <- read.xlsx("~/GoogleDrive/Biogeography/data/processed/counts_annotated.xlsx", 1)
annotation <- counts[,1:4]
librarysize <- colSums(counts[,5:21])
cpm <- data.frame(t((t(counts[,5:21]) / librarysize) * 10^6))
cpm$Identifier <- annotation$Identifier

asm.inner_asm.middle$Identifier <- rownames(asm.inner_asm.middle)
asm.inner_asm.middle <- merge(x=asm.inner_asm.middle, y=annotation,
                              by="Identifier")
asm.inner_asm.outer$Identifier <- rownames(asm.inner_asm.outer)
asm.inner_asm.outer <- merge(x=asm.inner_asm.outer, y=annotation, by = "Identifier")
min.inner_min.outer<-data.frame(min.inner_min.outer)
min.inner_min.outer$Identifier <- rownames(min.inner_min.outer)
min.inner_min.outer <- merge(x=min.inner_min.outer, y=annotation,
                             by="Identifier")

# filter
asm.inner.middle.upreg <- subset(asm.inner_asm.middle, FDR < 0.05 & 
                                   logFC > 1)
asm.inner.middle.downreg <- subset(asm.inner_asm.middle, FDR < 0.05 & 
                                   logFC < -1)


asm.inner.outer.upreg <- subset(asm.inner_asm.outer, FDR < 0.05 & 
                                   logFC > 1)
asm.inner.outer.downreg <- subset(asm.inner_asm.outer, FDR < 0.05 & 
                                     logFC < -1)

min.inner.outer.upreg <- subset(min.inner_min.outer, FDR < 0.05 & 
                                  logFC > 1)
min.inner.outer.downreg <- subset(min.inner_min.outer, FDR < 0.05 & 
                                  logFC < -1)



SHARED.OUTER.UPREG <- merge(x=min.inner.outer.upreg, asm.inner.outer.upreg,
                            by="Identifier", all=F)
SHARED.OUTER.DOWNREG<- merge(x=min.inner.outer.downreg, asm.inner.outer.downreg,
                            by="X", all=F)

ASM.SHARED.UPREG <- merge(x=asm.inner.middle.upreg, y=asm.inner.outer.upreg, by="Identifier")

ASM.SHARED.DOWNREG <- merge(x=asm.inner.middle.downreg, y=asm.inner.outer.downreg, by="Identifier")

'''
write.table(x=asm.inner.middle.downreg, file = "Downreg_ASM_Middle.txt", sep="\t", quote=FALSE)
write.table(x=asm.inner.outer.downreg, file = "Downreg_ASM_Outer.txt", sep="\t", quote=FALSE)
write.table(x=asm.inner.middle.upreg, file = "Upreg_ASM_Middle.txt", sep="\t", quote=FALSE)
write.table(x=asm.inner.outer.upreg, file = "Upreg_ASM_Outer.txt", sep="\t", quote=FALSE)
write.table(x=min.inner.outer.downreg, file = "Downreg_minA_Outer.txt", sep="\t", quote=FALSE)
write.table(x=min.inner.outer.upreg, file = "Upreg_minA_Outer.txt", sep="\t", quote=FALSE)
'''

ASM.SHARED.DOWNREG.counts <- merge(x=ASM.SHARED.DOWNREG,
                                   y=cpm,
                                   by.x="X", by.y="Identifier")

# plots of minA




