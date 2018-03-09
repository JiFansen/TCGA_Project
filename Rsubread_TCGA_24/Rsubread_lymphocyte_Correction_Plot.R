library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","LAML",
                      "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                      "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                      "UCS","UVM","ESCA","GBM","SKCM")
final.biospe.file <- data.frame()
for(i in 1:length(cancername.total)){
  biospecimen.name <- paste(cancername.total[i],".clean.no.zero")
  biospecimen.file <- read.table(file = biospecimen.name, header = T, sep = "\t")
  biospecimen.file <- biospecimen.file[8,]
  biospecimen.file <- t(biospecimen.file)
  colnames(biospecimen.file) <- "lymphocyte_infiltration"
  final.biospe.file <- rbind.data.frame(final.biospe.file,biospecimen.file)
}
final.biospe.file$sample.names <- rownames(final.biospe.file)

data.rsubread24 <- read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt",header=T,sep="\t",row.names=1)
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
label.biospe <- substr(as.character(rownames(final.biospe.file)),1,16)
label.rsubread <- substr(as.character(colnames(data.no.small)),1,16)
rsubread_2_biospe <- match(label.rsubread,label.biospe)
index.rsubread <- which(is.na(rsubread_2_biospe)==FALSE)
data.rsubread.subset <- data.no.small[,index.rsubread]
index.biospe <- rsubread_2_biospe[index.rsubread]
data.biospe.subset <- final.biospe.file[index.biospe,]
data.biospe.subset[,1] <- as.numeric(data.biospe.subset[,1])/100



xCellGene <- read.table(file="xCellGeneList",header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.no.small <- data.no.small[xcell.index,]