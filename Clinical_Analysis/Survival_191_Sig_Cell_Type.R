library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")

data.rsubread.correct <- read.table(file = 'data.subset',header = T, sep = '\t')
genename <- rownames(data.rsubread.correct)
expression.id <- colnames(data.rsubread.correct)

xCellScore <- read.table(file="xCellScore.subset",header=T,sep="\t")
antiGenScore <- read.table(file = 'antiGenScore.subset',header=T,sep="\t")

cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])
cancername <- unique(cancertype)

survival.patient.id <- read.table(file = "new+1.txt", header = F)
survival.patient.id <- as.character(survival.patient.id[,1])

common.id <- intersect(expression.id,survival.patient.id)
match.result <- match(expression.id, common.id)
big.index <- which(is.na(match.result)==FALSE)

data.subcancer <- data.rsubread.correct[,big.index]
subcancer.label <- cancertype[big.index]
subxCellScore <- xCellScore[,big.index]
subantiGenScore <- antiGenScore[,big.index]


cancershape.big <- c('A','B','C','D','E','F','G','H','I','K',"L",'M','N','O','P','R','S','T','U','V','W','X','Y','Z')
xcell.each.celltype.gene <- read.csv(file = "xCell_Type_Gene.csv", header = F, sep = ',', fill = T)
cellname <- unlist((strsplit(as.character(xcell.each.celltype.gene[,1]),"_")))
celllable <- c()
for(i in 1:length(cellname)){
  if(i%%3==1){
    celllable <- c(celllable,cellname[i])
  }
}
xcell.each.celltype.gene[,2] <- celllable

survival.sig.cells <- read.csv(file = "p_sig_score_name_new.csv",header = T, sep = ',')
survival.sig.cells <- as.character(survival.sig.cells[,2])

table.subcancer.lable <- as.data.frame(table(subcancer.label))

for(i in 1:(length(survival.sig.cells)-1)){
  cell.name <- survival.sig.cells[i]
  gene.subset <- xcell.each.celltype.gene[which(celllable==cell.name),3:202]
  gene.union <- character()
  for(j in 1:dim(gene.subset)[1]){
    gene.union <- union(gene.union,as.character(as.matrix(gene.subset[j,])))
  }
  gene.union <- intersect(gene.union,genename)
  genematch.result <- match(genename,gene.union)
  gene.index <- which(is.na(genematch.result)==FALSE)
  
  data.survival <- data.subcancer[gene.index,]
  
  data.survival.toal.SCALE <- t(scale(t(data.survival),center = TRUE, scale = TRUE))
  tsne <- Rtsne(t(data.survival.toal.SCALE), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
  data.pca <- t(data.survival.toal.SCALE)
  data.pr <- prcomp(data.pca)
  example_large_scale = SIMLR(X = data.survival.toal.SCALE, c = 12, cores.ratio = 0.5)
  
  g <- which(rownames(xCellScore)==cell.name)
  top <- quantile(as.numeric(subxCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[2]
  bottom <- quantile(as.numeric(subxCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[5]
  top.index <- which(as.numeric(subxCellScore[g,])<=top)
  bottom.index <- which(as.numeric(subxCellScore[g,])>bottom)
  between.index <- which(as.numeric(subxCellScore[g,])>top & as.numeric(subxCellScore[g,])<=bottom)
  
  m <- c()
  m[top.index] <- min(as.numeric(subxCellScore[g,]))
  m[bottom.index] <- max(as.numeric(subxCellScore[g,]))
  m[between.index] <- as.numeric(subxCellScore[g,])[between.index]
  
  dataset.tsne <- data.frame(tsne$Y,antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  rownames(dataset.tsne) <- colnames(subxCellScore)
  colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
  dataSet <- data.pr$x[,1:2]
  dataSet <- data.frame(dataSet, antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  dataset <- data.frame(example_large_scale$ydata, antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  rownames(dataset) <- colnames(subxCellScore)
  colnames(dataset)[c(1,2)] <- c("component1","component2")
  
  
  dataset.tsne.top <- dataset.tsne[top.index,]
  dataset.tsne.between <- dataset.tsne[between.index,]
  dataset.tsne.bottom <- dataset.tsne[bottom.index,]
  
  cancershape.small <- cancershape.big[1:length(unique(subcancer.label))]
  p.tsne <- ggplot()
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.tsne <- p.tsne+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  
  dataSet.top <- dataSet[top.index,]
  dataSet.between <- dataSet[between.index,]
  dataSet.bottom <- dataSet[bottom.index,]
  
  p.pca <- ggplot()
  p.pca <- p.pca+geom_point(data = dataSet.top, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.pca <- p.pca+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  p.pca <- p.pca+geom_point(data = dataSet.bottom, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'red')
  p.pca <- p.pca+geom_point(data = dataSet.between, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  dataset.top <- dataset[top.index,]
  dataset.between <- dataset[between.index,]
  dataset.bottom <- dataset[bottom.index,]
  
  p.simlr <- ggplot()
  p.simlr <- p.simlr+geom_point(data = dataset.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.simlr <- p.simlr+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
  p.simlr <- p.simlr+geom_point(data = dataset.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
  p.simlr <- p.simlr+geom_point(data = dataset.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  pdf(file = paste("data.subcancer.total.SCALE",cell.name,".pdf",sep = '') ,width=28, height=7)
  p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
  print(p)
  dev.off()
  
  data.survival.each.SCALE <- data.survival
  for(s in 1:dim(table.subcancer.lable)[1]){
    if(table.subcancer.lable[s,2]>2){
      cancer <- table.subcancer.lable[s,1]
      samples.Name <- which(subcancer.label==cancer)
      data.between <- data.survival.each.SCALE[,samples.Name]
      data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
      data.survival.each.SCALE[,samples.Name] <- data.between
    }
  }
  
  tsne <- Rtsne(t(data.survival.each.SCALE), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
  data.pca <- t(data.survival.each.SCALE)
  data.pr <- prcomp(data.pca)
  example_large_scale = SIMLR(X = data.survival.each.SCALE, c = 8, cores.ratio = 0.5)
  
  g <- which(rownames(xCellScore)==cell.name)
  top <- quantile(as.numeric(subxCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[2]
  bottom <- quantile(as.numeric(subxCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[5]
  top.index <- which(as.numeric(subxCellScore[g,])<=top)
  bottom.index <- which(as.numeric(subxCellScore[g,])>bottom)
  between.index <- which(as.numeric(subxCellScore[g,])>top & as.numeric(subxCellScore[g,])<=bottom)
  
  m <- c()
  m[top.index] <- min(as.numeric(subxCellScore[g,]))
  m[bottom.index] <- max(as.numeric(subxCellScore[g,]))
  m[between.index] <- as.numeric(subxCellScore[g,])[between.index]
  
  dataset.tsne <- data.frame(tsne$Y,antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  rownames(dataset.tsne) <- colnames(subxCellScore)
  colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
  dataSet <- data.pr$x[,1:2]
  dataSet <- data.frame(dataSet, antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  dataset <- data.frame(example_large_scale$ydata, antiGen.score = as.numeric(subantiGenScore), immune.score = m, cancertype = subcancer.label)
  rownames(dataset) <- colnames(subxCellScore)
  colnames(dataset)[c(1,2)] <- c("component1","component2")
  
  
  dataset.tsne.top <- dataset.tsne[top.index,]
  dataset.tsne.between <- dataset.tsne[between.index,]
  dataset.tsne.bottom <- dataset.tsne[bottom.index,]
  
  cancershape.small <- cancershape.big[1:length(unique(subcancer.label))]
  p.tsne <- ggplot()
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.tsne <- p.tsne+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  
  dataSet.top <- dataSet[top.index,]
  dataSet.between <- dataSet[between.index,]
  dataSet.bottom <- dataSet[bottom.index,]
  
  p.pca <- ggplot()
  p.pca <- p.pca+geom_point(data = dataSet.top, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.pca <- p.pca+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
  p.pca <- p.pca+geom_point(data = dataSet.bottom, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'red')
  p.pca <- p.pca+geom_point(data = dataSet.between, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  dataset.top <- dataset[top.index,]
  dataset.between <- dataset[between.index,]
  dataset.bottom <- dataset[bottom.index,]
  
  p.simlr <- ggplot()
  p.simlr <- p.simlr+geom_point(data = dataset.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
  p.simlr <- p.simlr+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape.small)
  p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
  p.simlr <- p.simlr+geom_point(data = dataset.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
  p.simlr <- p.simlr+geom_point(data = dataset.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  pdf(file = paste("data.subcancer.each.SCALE",cell.name,".pdf",sep = '') ,width=28, height=7)
  p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
  print(p)
  dev.off()
  
  
  
}


#i=length(survival.sig.cells)

         














































