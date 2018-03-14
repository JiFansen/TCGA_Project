library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])

#xCellGene <- c("SASH3","WAS","NCKAP1L","CD37","CD48","CD53")
xCellGene <- c("HLA-DRA","HLA-DPA1","HLA-DRB1","HLA-DPB1","CD74","HLA-DMA","HLA-DMB","HLA-DQA1")

#cancername <- c("OV","LUAD","LUSC","CESC","STAD","HNSC")
cancername <- c("LUAD","COAD","STAD","LAML","LGG")

#sub.cancertype.index <- which(cancertype=="OV" | cancertype=="LUAD" | cancertype=="LUSC" | cancertype=="CESC" | cancertype=="STAD" | cancertype=="HNSC")
sub.cancertype.index <- which(cancertype=="LUAD" | cancertype=="COAD" | cancertype=="STAD" | cancertype=="LAML" | cancertype=="LGG")

#cancershape <- c('A','B','C','D','E','F')
cancershape <- c('A','B','C','D','E')

#filename <- "xCell_6_gene_6_cancertype.pdf"
filename <- "antiGen_8_gene_5_cancertype.pdf"

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

data.xcell <- t(scale(t(data.xcell),center = TRUE, scale = TRUE))

#for(i in 1:length(cancername)){  #### Because GBM only has two samples.
#  cancer <- cancername[i]
#  samples.Name <- which(cancertype==cancer)
#  data.between <- data.xcell[,samples.Name]
#  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
#  data.xcell[,samples.Name] <- data.between
#}

#write.table(data.xcell, file = paste(filename,".txt", sep = ''), sep = '\t')

tsne <- Rtsne(t(data.xcell), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.xcell)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR(X = data.xcell, c = 8)

top <- quantile(as.numeric(xCellScore.subset[65,]), prob = c(0,0.2,0.4,0.6,0.8,1))[2]
bottom <- quantile(as.numeric(xCellScore.subset[65,]), prob = c(0,0.2,0.4,0.6,0.8,1))[5]
top.index <- which(as.numeric(xCellScore.subset[65,])<=top)
bottom.index <- which(as.numeric(xCellScore.subset[65,])>bottom)
between.index <- which(as.numeric(xCellScore.subset[65,])>top & as.numeric(xCellScore.subset[65,])<=bottom)

m <- c()
m[top.index] <- min(as.numeric(xCellScore.subset[65,]))
m[bottom.index] <- max(as.numeric(xCellScore.subset[65,]))
m[between.index] <- as.numeric(xCellScore.subset[65,])[between.index]


dataset.tsne <- data.frame(tsne$Y,antiGen.score = antiGenScore.subset, immune.score1 = as.numeric(xCellScore.subset[65,]), immune.score = m, cancertype = cancertype)
rownames(dataset.tsne) <- colnames(xCellScore.subset)
colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet, antiGen.score = antiGenScore.subset, immune.score2 = as.numeric(xCellScore.subset[65,]), immune.score = m, cancertype = cancertype)
dataset <- data.frame(example_large_scale$ydata, antiGen.score = antiGenScore.subset, immune.score3 = as.numeric(xCellScore.subset[65,]), immune.score = m, cancertype = cancertype)
rownames(dataset) <- colnames(xCellScore.subset)
colnames(dataset)[c(1,2)] <- c("component1","component2")


dataset.tsne.top <- dataset.tsne[top.index,]
dataset.tsne.between <- dataset.tsne[between.index,]
dataset.tsne.bottom <- dataset.tsne[bottom.index,]

p.tsne <- ggplot()
p.tsne <- p.tsne+geom_point(data = dataset.tsne.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
p.tsne <- p.tsne+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
p.tsne <- p.tsne+geom_point(data = dataset.tsne.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
p.tsne <- p.tsne+geom_point(data = dataset.tsne.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')

dataSet.top <- dataSet[top.index,]
dataSet.between <- dataSet[between.index,]
dataSet.bottom <- dataSet[bottom.index,]

p.pca <- ggplot()
p.pca <- p.pca+geom_point(data = dataSet.top, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'blue')
p.pca <- p.pca+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
p.pca <- p.pca+geom_point(data = dataSet.bottom, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype), col = 'red')
p.pca <- p.pca+geom_point(data = dataSet.between, aes(x=PC1,y=PC2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')

dataset.top <- dataset[top.index,]
dataset.between <- dataset[between.index,]
dataset.bottom <- dataset[bottom.index,]

p.simlr <- ggplot()
p.simlr <- p.simlr+geom_point(data = dataset.top, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'blue')
p.simlr <- p.simlr+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))+theme(legend.position="none")
p.simlr <- p.simlr+geom_point(data = dataset.bottom, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype), col = 'red')
p.simlr <- p.simlr+geom_point(data = dataset.between, aes(x=component1,y=component2,size = antiGen.score, shape = cancertype, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')


pdf(file=paste(filename, "_ALL_SCALE", sep = '') ,width=28, height=7)
p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
print(p)
dev.off()

write.table(dataset.tsne, file = paste(filename,"-tsne", sep = ''), sep = '\t')
write.table(dataSet, file = paste(filename,"-pca", sep = ''), sep = '\t')
write.table(dataset, file = paste(filename,"-simlr", sep = ''), sep = '\t')
write.table(cancertype, file = paste(filename,"cancertype", sep = ''), sep = '\t')


dataset <- read.table(file = 'antiGen_8_gene_5_cancertype.pdf-simlr',header = T,sep = '\t')
top <- quantile(as.numeric(dataset$immune.score3), prob = c(0,0.2,0.4,0.6,0.8,1))[2]
bottom <- quantile(as.numeric(dataset$immune.score3), prob = c(0,0.2,0.4,0.6,0.8,1))[5]
top.index <- which(as.numeric(dataset$immune.score3)<=top)
bottom.index <- which(as.numeric(dataset$immune.score3)>bottom)
between.index <- which(as.numeric(dataset$immune.score3)>top & as.numeric(dataset$immune.score3)<=bottom)



p.simlr+geom_vline(aes(xintercept=-45))+geom_vline(aes(xintercept=-32))+geom_hline(aes(yintercept=28))+geom_hline(aes(yintercept=37))
data.1 <- dataset[which(dataset$component1>=-45 & dataset$component1<=-32),]
data.2 <- data.1[which(data.1$component2>=28 & data.1$component2<=37),]
Name <- rownames(data.2)
write.table(Name,"antiGen_8_gene_5_cancertype_Immune_Top.txt", sep = '\t')

p.simlr+geom_vline(aes(xintercept=8))+geom_vline(aes(xintercept=25))+geom_hline(aes(yintercept=-50))+geom_hline(aes(yintercept=-23))
data.1 <- dataset[which(dataset$component1>=8 & dataset$component1<=25),]
data.2 <- data.1[which(data.1$component2>=-50 & data.1$component2<=-23),]
Name <- rownames(data.2)
write.table(Name,"antiGen_8_gene_5_cancertype_Immune_bottom.txt", sep = '\t')
