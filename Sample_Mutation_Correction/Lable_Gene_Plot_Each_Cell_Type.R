library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")

featureName <- "xCellGeneList"
#featureName <- "AntigenProcessing_KEGG.txt"

data.rsubread.correct <- read.table(file = 'data.subset',header = T, sep = '\t')
sampleName <- colnames(data.rsubread.correct)
genename <- rownames(data.rsubread.correct)


cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])

antiGenScore <- read.table(file = 'antiGenScore.subset',header=T,sep="\t")
SASH3Score <- as.numeric(data.rsubread.correct[which(genename=="SASH3"),])

xCellGene <- read.table(file=featureName,header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)

cancer.total <- as.character(unique(cancertype))
for(n in 1:length(cancer.total)){
  inputName <- cancer.total[n]
  brca.index <- which(cancertype==inputName)
  data.brca <- data.rsubread.correct[,brca.index]
  antiGenScore.brca <- antiGenScore[,brca.index]
  SASH3Score.brca <- SASH3Score[brca.index]
  data.brca <- data.brca[xcell.index,]
  data.brca <- t(scale(t(data.brca)))
  
  tsne <- Rtsne(t(data.brca), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
  data.pca <- t(data.brca)
  data.pr <- prcomp(data.pca)
  example_large_scale = SIMLR(X = data.brca, c = 8, cores.ratio = 0.5)
  
  dataset.tsne <- data.frame(tsne$Y,antiGen.score = as.numeric(antiGenScore.brca[1,]), SASH3.score = SASH3Score.brca)
  colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
  dataSet <- data.pr$x[,1:2]
  dataSet <- data.frame(dataSet, antiGen.score = as.numeric(antiGenScore.brca[1,]), SASH3.score = SASH3Score.brca)
  dataset <- data.frame(example_large_scale$ydata, antiGen.score = as.numeric(antiGenScore.brca[1,]), SASH3.score = SASH3Score.brca)
  colnames(dataset)[c(1,2)] <- c("component1","component2")
  
  p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,size = SASH3.score, col = antiGen.score))
  p.tsne <- p.tsne+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))
  p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient2(low = 'blue', midpoint = median(dataset.tsne$antiGen.score), mid = 'green', high = 'red')
  
  
  p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,size = SASH3.score, col = antiGen.score))
  p.pca <- p.pca+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  p.pca <- p.pca+scale_colour_gradient2(low = 'blue', midpoint = median(dataSet$antiGen.score), mid = 'green', high = 'red')
  
  
  p.simlr <- ggplot(dataset,aes(x=component1,y=component2,size = SASH3.score, col = antiGen.score))
  p.simlr <- p.simlr+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))
  p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient2(low = 'blue', midpoint = median(dataset$antiGen.score), mid = 'green', high = 'red')
  
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/SASH3_Lable_Each_Cancer_Type")
  pdf(file=paste(inputName,"_",featureName,"_Label_SASH3.pdf",sep=''),width=28, height=7)
  p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
  print(p)
  dev.off()
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/")
}