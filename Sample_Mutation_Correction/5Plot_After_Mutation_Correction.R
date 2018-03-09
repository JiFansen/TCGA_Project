data.final <- read.table(file = "data.subset", sep = "\t", header = T)
genename <- rownames(data.final)
library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
xCellGene <- read.table(file="AntigenProcessing_KEGG.txt",header=T,sep="\t")

cancertype <- read.table(file="CancerLable",header=T,sep="\t")
cancername <- unique(as.character(cancertype[,1]))
#cancertype$sampleBarcode <- substr(as.character(cancertype$V1),1,16)

xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.final <- data.final[xcell.index,]
for(i in 1:length(cancername)){
  cancer <- cancername[i]
  samples.Name <- which(as.character(cancertype[,1])==cancer)
  data.between <- data.final[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.final[,samples.Name] <- data.between
}
xCellScore <- read.table(file="xCellScore.subset",header=T,sep="\t")
tsne <- Rtsne(t(data.final), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.final)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR_Large_Scale(X = data.final, c = 8, kk = 10)

dataset.tsne <- data.frame(tsne$Y,cancertype=as.character(cancertype[,1]),score = as.numeric(xCellScore[65,]))
colnames(dataset.tsne) <- c("component1","component2","cancertype","score")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet,cancertype=cancertype[,1],score = as.numeric(xCellScore[65,]))
dataset <- data.frame(example_large_scale$ydata,cancertype=as.character(cancertype[,1]),score = as.numeric(xCellScore[65,]))
colnames(dataset) <- c("component1","component2","cancertype","score")
p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,col=cancertype,size=score))
p.tsne <- p.tsne+geom_point()+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))
p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")
p.simlr <- ggplot(dataset,aes(x=component1,y=component2,col=cancertype,size=score))
p.simlr <- p.simlr+geom_point()+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))
p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")
p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,col=cancertype,size=score))
p.pca <- p.pca+geom_point()+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
#picture.name <- paste("xCellScore_",rownames(xCellScore)[i],".png",sep="")
png(file="Mutation_Correction_9624.png",width=200, height=200, units='mm',res=1600)
ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
dev.off()
