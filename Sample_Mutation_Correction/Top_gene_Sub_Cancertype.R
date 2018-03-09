library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])

xCellGene <- c("SASH3","WAS","NCKAP1L","CD37","CD48","CD53")
#xCellGene <- c("HLA-DRA","HLA-DPA1","HLA-DRB1","HLA-DPB1","CD74","HLA-DMA","HLA-DMB","HLA-DQA1")

cancername <- c("OV","LUAD","LUSC","CESC","STAD","HNSC")
#cancername <- c("LUAD","COAD","STAD","LAML","LGG")

sub.cancertype.index <- which(cancertype=="OV" | cancertype=="LUAD" | cancertype=="LUSC" | cancertype=="CESC" | cancertype=="STAD" | cancertype=="HNSC")
#sub.cancertype.index <- which(cancertype=="LUAD" | cancertype=="COAD" | cancertype=="STAD" | cancertype=="LAML" | cancertype=="LGG")

cancershape <- c('A','B','C','D','E','F')
#cancershape <- c('A','B','C','D','E')

filename <- "xCell_6_gene_6_cancertype.pdf"
#filename <- "antiGen_8_gene_5_cancertype.pdf"

data.rsubread.correct <- read.table(file = 'data.subset',header = T, sep = '\t')
sampleName <- colnames(data.rsubread.correct)
genename <- rownames(data.rsubread.correct)


xCellScore <- read.table(file="xCellScore.subset",header=T,sep="\t")
antiGenScore <- read.table(file = 'antiGenScore.subset',header=T,sep="\t")

xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.xcell <- data.rsubread.correct[xcell.index,]

data.xcell <- data.xcell[,sub.cancertype.index]
cancertype <- cancertype[sub.cancertype.index]
xCellScore.subset <- xCellScore[,sub.cancertype.index]
antiGenScore.subset <- as.numeric(antiGenScore[1,sub.cancertype.index])

for(i in 1:length(cancername)){  #### Because GBM only has two samples.
  cancer <- cancername[i]
  samples.Name <- which(cancertype==cancer)
  data.between <- data.xcell[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.xcell[,samples.Name] <- data.between
}


tsne <- Rtsne(t(data.xcell), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.xcell)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR(X = data.xcell, c = 8)

dataset.tsne <- data.frame(tsne$Y,antiGen.score = antiGenScore.subset, immune.score = as.numeric(xCellScore.subset[65,]))
colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet, antiGen.score = antiGenScore.subset, immune.score = as.numeric(xCellScore.subset[65,]))
dataset <- data.frame(example_large_scale$ydata, antiGen.score = antiGenScore.subset, immune.score = as.numeric(xCellScore.subset[65,]))
colnames(dataset)[c(1,2)] <- c("component1","component2")

p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,size = antiGen.score, col = immune.score, shape = cancertype))
p.tsne <- p.tsne+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+scale_shape_manual(values=cancershape)
p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient2(low = 'blue', midpoint = median(dataset.tsne$immune.score), mid = 'green', high = 'red')


p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,size = antiGen.score, col = immune.score, shape = cancertype))
p.pca <- p.pca+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+scale_shape_manual(values=cancershape)
p.pca <- p.pca+scale_colour_gradient2(low = 'blue', midpoint = median(dataSet$immune.score), mid = 'green', high = 'red')


p.simlr <- ggplot(dataset,aes(x=component1,y=component2,size = antiGen.score, col = immune.score, shape = cancertype))
p.simlr <- p.simlr+geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+scale_shape_manual(values=cancershape)
p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient2(low = 'blue', midpoint = median(dataset$immune.score), mid = 'green', high = 'red')

pdf(file=filename,width=28, height=7)
p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
print(p)
dev.off()
