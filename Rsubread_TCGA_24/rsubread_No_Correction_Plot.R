library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
data.rsubread24 <- read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt",header=T,sep="\t",row.names=1)
############################################## xCell Genes or any other genelists. ################################################################
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
#xCellGene <- read.table(file="xCellGeneList",header=T,sep="\t")
xCellGene <- read.table(file="AntigenProcessing_KEGG.txt",header=T,sep="\t")
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
############################################################################################################################
xCellScore <- read.table(file="xCellScore.txt",header=T,sep="\t")

tsne <- Rtsne(t(data.no.small), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.no.small)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)

setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/")
for(i in 1:dim(xCellScore)[1]){
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/point_cell_type_score_antigen")
  dataset.tsne <- data.frame(tsne$Y,cancertype=as.character(cancertype[,2]),score = as.numeric(xCellScore[i,]))
  colnames(dataset.tsne) <- c("component1","component2","cancertype","score")
  dataSet <- data.pr$x[,1:2]
  dataSet <- data.frame(dataSet,cancertype=cancertype[,2],score = as.numeric(xCellScore[i,]))
  dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]),score = as.numeric(xCellScore[i,]))
  colnames(dataset) <- c("component1","component2","cancertype","score")
  p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,col=cancertype,size=score))
  p.tsne <- p.tsne+geom_point()+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))
  p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")
  p.simlr <- ggplot(dataset,aes(x=component1,y=component2,col=cancertype,size=score))
  p.simlr <- p.simlr+geom_point()+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))
  p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")
  p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,col=cancertype,size=score))
  p.pca <- p.pca+geom_point()+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  picture.name <- paste("xCellScore_",rownames(xCellScore)[i],".png",sep="")
  png(file=picture.name,width=200, height=200, units='mm',res=1600)
  ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
  dev.off()
}
##########################  Label the samples using each of the 67 scores. #####################################
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/")
example_large_scale = SIMLR_Large_Scale(X = data.no.small, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]))
colnames(dataset) <- c("component1","component2","cancertype")
for(i in 1:dim(xCellScore)[1]){
  picture.name <- paste("xCellScore_",rownames(xCellScore)[i],".png",sep="")
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/point_cell_type_score")
  png(file=picture.name,width=200, height=200, units='mm',res=1600)
  size <- as.numeric(xCellScore[i,])
  p <- ggplot(dataset,aes(x=component1,y=component2,col=cancertype,size=size))
  p <- p+geom_point()+labs(title=picture.name)+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  print(p)
  dev.off()
}

################################ Using the Genelists of each cell types to cluster the samples. ######################
################################ GeneList Preparing. #################################################################
setwd("C:/Users/lenovo/Desktop")
total.gene <- read.csv(file="xCell_CellType.csv",header=F,fill=T)
celltype <- character()
for(i in 1:dim(total.gene)[1]){
  celltype[i] <- unlist(strsplit(as.character(total.gene[,1]),"_")[i])[1]
}
total.gene$celltype <- celltype
total.gene <- total.gene[,-c(1,2)]
uniq.cell.type <- unique(as.character(total.gene$celltype))
for (i in 1:length(uniq.cell.type)){
  data.between <- total.gene[total.gene$celltype==uniq.cell.type[i],]
  data.between <- data.between[,-ncol(data.between)]
  gene.each.cell.type <- character()
  for(j in 1:dim(data.between)[1]){
    gene.each.cell.type <- union(gene.each.cell.type,as.character(as.matrix(data.between[j,])))
  }
  nothing.index <- which(gene.each.cell.type=="")
  gene.each.cell.type <- gene.each.cell.type[-nothing.index]
  setwd("C:/Users/lenovo/Desktop/Gene_Each_Cell_Type/")
  write.table(gene.each.cell.type,paste(uniq.cell.type[i],"_GeneList",sep=""),sep="\t")
}

############################### Using the Genelists of each cell types to cluster the samples. #######################
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
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Gene_Each_Cell_Type")
uniq.cell.type <- read.table(file = "uniq.cell.type",sep="\t")
uniq.cell.type <- as.character(uniq.cell.type[,1])
for(m in 1:length(uniq.cell.type)){
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Gene_Each_Cell_Type")
  xCellGene <- read.table(file=paste(uniq.cell.type[m],"_GeneList",sep=""),header=T,sep="\t")
  xCellGene <- as.character(xCellGene[,1])
  xcellGeneList <- intersect(xCellGene,genename)
  xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
  data.subset <- data.no.small[xcell.index,]
  for(i in 1:length(cancername)){
    cancer <- cancername[i]
    samples.Name <- which(as.character(cancertype[,2])==cancer)
    data.between <- data.subset[,samples.Name]
    data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
    data.subset[,samples.Name] <- data.between
  }
  example_large_scale = SIMLR_Large_Scale(X = data.subset, c = 8, kk = 10)
  dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,2]))
  colnames(dataset) <- c("component1","component2","cancertype")
  filename <- paste(uniq.cell.type[m],".png",sep="")
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Plot_Each_Cell_Type")
  png(file=filename,width=200, height=200, units='mm',res=1600)
  p <- ggplot(dataset,aes(x=component1,y=component2,col=cancertype))
  p <- p+geom_point()+labs(title=filename)+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  print(p)
  dev.off()
}

######################################## Jing Zhe's test. xCell Genes. #########################################################
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
