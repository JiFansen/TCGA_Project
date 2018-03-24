setwd("/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/correct")
i = 1
for (fileName in list.files("/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/correct", "*correct")){
  data.between <- read.table(file = fileName, header = T, sep = '\t')
  if(i == 1){
    data.result <- data.between
  }else{
    data.result <- rbind.data.frame(data.result, data.between)
  }
  i = i+1
}
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
write.table(data.result,"total.score.table", sep = '\t')

function.score <- read.table(file = "total.score.table", header = T, sep = '\t')
sampleName <- read.table(file = '/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/correct/sample', header = T, sep = '\t')
sampleName <- as.character(sampleName[,1])
sampleName <- gsub("\\.","-",sampleName)
sampleName <- substr(sampleName,1,16)

geneName <- read.table(file = '/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/correct/Name', header = T, sep = '\t')
geneName <- as.character(geneName[,1])

function.matrix <- matrix(1, nrow = length(geneName), ncol = length(sampleName))
rownames(function.matrix) <- geneName
colnames(function.matrix) <- sampleName

match.result <- match(as.character(function.score$sample), colnames(function.matrix))
Sample.index <- which(is.na(match.result)==FALSE)

new.function.score <- function.score[Sample.index,]

for(i in 1:dim(new.function.score)[1]){
  gene.index <- which(rownames(function.matrix)==as.character(new.function.score[i,1]))
  sample.index <- which(colnames(function.matrix)==as.character(new.function.score[i,2]))
  function.matrix[gene.index,sample.index] <- as.numeric(new.function.score[i,3])
}

write.table(function.matrix,"function.matrix", sep = '\t')
###############################################################################################################
library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")
library("metap")



# Feature gene List.
featureName <- "AntigenProcessing_GeneList"
cell.name <- unlist(strsplit(featureName,"_"))[1]
# Read the Score data matrix.
xCellScore <- read.table(file="AntigenProcessing.subset",header=T,sep="\t")
# Read the corrected expression data matrix.
function.matrix <- read.table(file = "function.matrix", header = T, sep = '\t')
sampleName <- colnames(function.matrix)
geneName <- rownames(function.matrix)

# Read the feature gene name.
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Signiture_Gene")
xCellGene <- read.table(file=featureName,header=T,sep="\t")
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)

# Read the cancer label matrix.
cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])
cancername <- unique(cancertype)

cancershape <- c('A','B','C','D','E','F','G','H','I','K',"L",'M','N','O','P','R','S','T','U','V','W','X','Y','Z')
for(i in 1:length(cancername)){
  inputName <- cancername[i]
  brca.index <- which(cancertype==inputName)
  data.brca <- function.matrix[,brca.index]
  subxCellScore <- xCellScore[,brca.index]
  data.brca <- data.brca[xcell.index,]
  data.brca <- t(scale(t(data.brca)))
  
  tsne <- Rtsne(t(data.brca), dims = 2, perplexity=5, verbose=TRUE, max_iter = 500)
  data.pca <- t(data.brca)
  data.pr <- prcomp(data.pca)
  example_large_scale = SIMLR(X = data.brca, c = 5, cores.ratio = 0.5)
  ranks = SIMLR_Feature_Ranking(A=example_large_scale$S,X = data.brca)
  genes.rank <- data.frame(gene = rownames(data.brca)[ranks$aggR], p.value = as.numeric(ranks$pval))
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Each_Cell_Type_Gene_Rank")
  write.table(genes.rank, paste("functionScore",inputName,cell.name,"_genes.rank",sep = ''),sep = '\t')
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
  if(i==1){
    GeneRank.matrix <- genes.rank[order(genes.rank$gene),]
  }else{
    GeneRank.matrix <- cbind.data.frame(GeneRank.matrix, genes.rank[order(genes.rank$gene),])
  }
  
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
  
  
  dataset.tsne <- data.frame(tsne$Y, immune.score = m)
  rownames(dataset.tsne) <- colnames(subxCellScore)
  colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
  dataSet <- data.pr$x[,1:2]
  dataSet <- data.frame(dataSet, immune.score = m)
  dataset <- data.frame(example_large_scale$ydata, immune.score = m)
  rownames(dataset) <- colnames(subxCellScore)
  colnames(dataset)[c(1,2)] <- c("component1","component2")
  
  dataset.tsne.top <- dataset.tsne[top.index,]
  dataset.tsne.between <- dataset.tsne[between.index,]
  dataset.tsne.bottom <- dataset.tsne[bottom.index,]
  
  p.tsne <- ggplot()
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.top, aes(x=component1,y=component2), col = 'blue')
  p.tsne <- p.tsne+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+theme(legend.position="none")
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.bottom, aes(x=component1,y=component2), col = 'red')
  p.tsne <- p.tsne+geom_point(data = dataset.tsne.between, aes(x=component1,y=component2, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  dataSet.top <- dataSet[top.index,]
  dataSet.between <- dataSet[between.index,]
  dataSet.bottom <- dataSet[bottom.index,]
  
  p.pca <- ggplot()
  p.pca <- p.pca+geom_point(data = dataSet.top, aes(x=PC1,y=PC2), col = 'blue')
  p.pca <- p.pca+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))
  p.pca <- p.pca+geom_point(data = dataSet.bottom, aes(x=PC1,y=PC2), col = 'red')
  p.pca <- p.pca+geom_point(data = dataSet.between, aes(x=PC1,y=PC2,col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  dataset.top <- dataset[top.index,]
  dataset.between <- dataset[between.index,]
  dataset.bottom <- dataset[bottom.index,]
  
  p.simlr <- ggplot()
  p.simlr <- p.simlr+geom_point(data = dataset.top, aes(x=component1,y=component2), col = 'blue')
  p.simlr <- p.simlr+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
  p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+theme(legend.position="none")
  p.simlr <- p.simlr+geom_point(data = dataset.bottom, aes(x=component1,y=component2), col = 'red')
  p.simlr <- p.simlr+geom_point(data = dataset.between, aes(x=component1,y=component2, col = immune.score))+scale_color_gradient(low = 'blue', high = 'red')
  
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Signiture_Gene_Plot")
  pdf(file = paste("functionScore",inputName,"_",cell.name,".pdf",sep = '') ,width=28, height=7)
  p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
  print(p)
  dev.off()
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
}
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Each_Cell_Type_Gene_Rank")
write.table(GeneRank.matrix, file = paste("functionScore",cell.name,"_Total_GeneRank.matrix",sep = ''), sep = '\t')

meta.p <- as.numeric(apply(GeneRank.matrix[,seq(2,48,2)],1,FUN=function(x){logitp(x)$p}))
logp.FDR<-p.adjust(meta.p,method="fdr",n=length(meta.p))
#meta.p <- as.numeric(1-pchisq(apply(log(GeneRank.matrix[,seq(2,48,2)]),1,sum)*(-2),2*length(cancername)))
genes <- as.character(GeneRank.matrix[,1])
meta.result <- data.frame(Gene = genes, p.value = meta.p,FDR = logp.FDR)
write.table(meta.result, file = paste("functionScore",cell.name,"meta.p.value"), sep = '\t')
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")

data.all.sample <- function.matrix
for(i in 1:length(cancername)){  
  cancer <- cancername[i]
  samples.Name <- which(cancertype==cancer)
  data.between <- data.xcell[,samples.Name]
  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
  data.all.sample[,samples.Name] <- data.between
}

tsne <- Rtsne(t(data.all.sample), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.all.sample)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR_Large_Scale(X = data.brca, c = 8, kk=10)

g <- which(rownames(xCellScore)==cell.name)
top <- quantile(as.numeric(xCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[2]
bottom <- quantile(as.numeric(xCellScore[g,]), prob = c(0,0.2,0.4,0.6,0.8,1))[5]
top.index <- which(as.numeric(xCellScore[g,])<=top)
bottom.index <- which(as.numeric(xCellScore[g,])>bottom)
between.index <- which(as.numeric(xCellScore[g,])>top & as.numeric(xCellScore[g,])<=bottom)
m <- c()
m[top.index] <- min(as.numeric(xCellScore[g,]))
m[bottom.index] <- max(as.numeric(xCellScore[g,]))
m[between.index] <- as.numeric(xCellScore[g,])[between.index]

dataset.tsne <- data.frame(tsne$Y, immune.score = m, cancertype = cancertype)
rownames(dataset.tsne) <- colnames(xCellScore)
colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet, immune.score = m, cancertype = cancertype)
dataset <- data.frame(example_large_scale$ydata, immune.score = m, cancertype = cancertype)
rownames(dataset) <- colnames(xCellScore)
colnames(dataset)[c(1,2)] <- c("component1","component2")

dataset.tsne.top <- dataset.tsne[top.index,]
dataset.tsne.between <- dataset.tsne[between.index,]
dataset.tsne.bottom <- dataset.tsne[bottom.index,]

p.tsne <- ggplot()
p.tsne <- p.tsne+geom_point(data = dataset.tsne.top, aes(x=component1,y=component2, shape = cancertype), col = 'blue')
p.tsne <- p.tsne+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.tsne <- p.tsne+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))+theme(legend.position="none")
p.tsne <- p.tsne+geom_point(data = dataset.tsne.bottom, aes(x=component1,y=component2, shape = cancertype), col = 'red')
p.tsne <- p.tsne+geom_point(data = dataset.tsne.between, aes(x=component1,y=component2, col = immune.score, shape = cancertype))+scale_color_gradient(low = 'blue', high = 'red')

dataSet.top <- dataSet[top.index,]
dataSet.between <- dataSet[between.index,]
dataSet.bottom <- dataSet[bottom.index,]

p.pca <- ggplot()
p.pca <- p.pca+geom_point(data = dataSet.top, aes(x=PC1,y=PC2, shape = cancertype), col = 'blue')
p.pca <- p.pca+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.pca <- p.pca+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))
p.pca <- p.pca+geom_point(data = dataSet.bottom, aes(x=PC1,y=PC2, shape = cancertype), col = 'red')
p.pca <- p.pca+geom_point(data = dataSet.between, aes(x=PC1,y=PC2,col = immune.score, shape = cancertype))+scale_color_gradient(low = 'blue', high = 'red')

dataset.top <- dataset[top.index,]
dataset.between <- dataset[between.index,]
dataset.bottom <- dataset[bottom.index,]

p.simlr <- ggplot()
p.simlr <- p.simlr+geom_point(data = dataset.top, aes(x=component1,y=component2, shape = cancertype), col = 'blue')
p.simlr <- p.simlr+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+scale_shape_manual(values=cancershape)
p.simlr <- p.simlr+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))+theme(legend.position="none")
p.simlr <- p.simlr+geom_point(data = dataset.bottom, aes(x=component1,y=component2, shape = cancertype), col = 'red')
p.simlr <- p.simlr+geom_point(data = dataset.between, aes(x=component1,y=component2, col = immune.score, shape = cancertype))+scale_color_gradient(low = 'blue', high = 'red')

setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Signiture_Gene_Plot")
pdf(file = paste("functionScore_Total_CancerType","_",cell.name,".pdf",sep = '') ,width=28, height=7)
p <- ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
print(p)
dev.off()
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")

