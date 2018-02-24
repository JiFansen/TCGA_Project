library('preprocessCore')
library(SIMLR)
library(igraph)
library(ggplot2)
##################################### Using All the Genes. ###########################################################
data.rsubread24 <- read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt",header=T,sep="\t",rownames=1)
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
cancertype <- read.table(file="cancertype",header=F,sep="\t")
cancername <- unique(as.character(cancertype[,2]))
for(i in 1:length(cancername)){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype[,2])==cancer)
  data.between <- data.no.small[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.no.small[,samples.Name] <- data.between
}
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]))
colnames(dataset) <- c("component1","component2","cancertype")
png(file="Rsubread_No_Correction_Plot_All_Gene.png",width=300, height=100, units='mm',res=1500)
ggplot(dataset,aes(x=component1,y=component2,col=cancertype))+geom_point()+labs(title="Rsubread_No_Correction_Plot_All_Gene.png")+theme(plot.title=element_text(hjust=0.5))
dev.off()

############################## Using the xCell Genelists. ###########################################################
############################## After xcellGene extraction, then qq normalization.####################################
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
xCellGene <- read.table(file="xCellGeneList",header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)

data.no.small <- data.no.small[xcell.index,]
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- xcellGeneList
colnames(data.no.small) <- sampleName
cancertype <- read.table(file="cancertype",header=F,sep="\t")
cancername <- unique(as.character(cancertype[,2]))
for(i in 1:length(cancername)){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype[,2])==cancer)
  data.between <- data.no.small[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.no.small[,samples.Name] <- data.between
}
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]))
colnames(dataset) <- c("component1","component2","cancertype")
png(file="Rsubread_No_Correction_Plot_xCell_Gene_1.png",width=300, height=100, units='mm',res=1500)
ggplot(dataset,aes(x=component1,y=component2,col=cancertype))+geom_point()+labs(title="Rsubread_No_Correction_Plot_xCell_Gene_1.png")+theme(plot.title=element_text(hjust=0.5))
dev.off()

############################## Using the xCell Genelists. ###########################################################
############################## First, qq normalization, then xcellGene extraction.####################################
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
xCellGene <- read.table(file="xCellGeneList",header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
data.no.small <- data.no.small[xcell.index,]
cancertype <- read.table(file="cancertype",header=F,sep="\t")
cancername <- unique(as.character(cancertype[,2]))
for(i in 1:length(cancername)){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype[,2])==cancer)
  data.between <- data.no.small[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.no.small[,samples.Name] <- data.between
}
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]))
colnames(dataset) <- c("component1","component2","cancertype")
png(file="Rsubread_No_Correction_Plot_xCell_Gene_1.png",width=300, height=100, units='mm',res=1500)
ggplot(dataset,aes(x=component1,y=component2,col=cancertype))+geom_point()+labs(title="Rsubread_No_Correction_Plot_xCell_Gene_1.png")+theme(plot.title=element_text(hjust=0.5))
dev.off()

############################################### Jing Zhe's test. All Genes. #####################################################
library('preprocessCore')
library(SIMLR)
library(igraph)
library(ggplot2)
data.715 <- read.table(file="tumor_715.txt",header=T,sep="\t",row.names=1)
data.no.small <- data.715[which(rowSums(data.715<1)!=dim(data.715)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
cancertype <- substr(colnames(data.715),1,4)
cancername <- unique(as.character(cancertype))
for(i in 1:16){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype)==cancer)
  data.between <- data.no.small[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.no.small[,samples.Name] <- data.between
}
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype))
colnames(dataset) <- c("component1","component2","cancertype")
png(file="Tumor_715_All_Gene.png",width=300, height=100, units='mm',res=1500)
ggplot(dataset,aes(x=component1,y=component2,col=cancertype))+geom_point()+labs(title="Tumor_715_All_Gene")+theme(plot.title=element_text(hjust=0.5))
dev.off()
######################################## Jing Zhe's test. xCell Genes. #########################################################
data.no.small <- data.715[which(rowSums(data.715<1)!=dim(data.715)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
cancertype <- substr(colnames(data.715),1,4)
cancername <- unique(as.character(cancertype))
xCellGene <- read.table(file="xCellGeneList",header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.no.small <- data.no.small[xcell.index,]
for(i in 1:16){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype)==cancer)
  data.between <- data.no.small[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.no.small[,samples.Name] <- data.between
}
xCellScore <- read.table(file="xCellScore",header=T,sep="\t")
size <- as.numeric(xCellScore[65,])
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype))
colnames(dataset) <- c("component1","component2","cancertype")
png(file="Tumor_715_xCell_Gene.png",width=200, height=200, units='mm',res=1600)
p <- ggplot(dataset,aes(x=component1,y=component2,col=cancertype,size=size))
p+geom_point()+labs(title="Tumor_715_xCell_Gene")+theme(plot.title=element_text(hjust=0.5))
dev.off()
