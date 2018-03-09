library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")


inputName <- 'UCEC'
featureName <- "xCellGeneList"
#featureName <- "AntigenProcessing_KEGG.txt"

data.rsubread <- read.table(file = 'GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt',header = T, sep = '\t', row.names = 1)
data.no.small <- data.rsubread[which(rowSums(data.rsubread<1)!=dim(data.rsubread)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName


cancertype <- read.table(file = 'cancertype',header = T, sep = '\t')
cancertype <- as.character(cancertype[,2])



brca.index <- which(cancertype==inputName)

data.brca <- data.no.small[,brca.index]
write.table(data.brca,paste(inputName,"_rsubread_No_Mutation_Correct.txt",sep=), sep = '\t')


xCellScore <- read.table(file="xCellScore.txt",header=T,sep="\t")
xCellScore.brca <- xCellScore[,brca.index]

xCellGene <- read.table(file=featureName,header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.brca <- data.brca[xcell.index,]


biospecimen.name <- paste(inputName,".clean.no.zero")
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Slides_Data")
biospecimen.file <- read.table(file = biospecimen.name, header = T, sep = "\t")
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
biospecimen.file <- biospecimen.file[8,]

expression.sample.name <- substr(colnames(data.brca),1,16)
biospeci.sample.name <- substr(colnames(biospecimen.file),1,16)

match.vector <- match(expression.sample.name,biospeci.sample.name)
expression.index <- which(is.na(match.vector)==FALSE)
biospeci.index <- match.vector[expression.index]

data.brca <- data.brca[,expression.index]
biospecimen.file <- biospecimen.file[,biospeci.index]
xCellScore.brca <- xCellScore.brca[,expression.index]

data.brca <- t(scale(t(data.brca)))

tsne <- Rtsne(t(data.brca), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.brca)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR(X = data.brca, c = 8, cores.ratio = 0.5)


ranks = SIMLR_Feature_Ranking(A=example_large_scale$S,X = data.brca)
genes.rank <- data.frame(gene = rownames(data.brca)[ranks$aggR], p.value = ranks$pval)
write.table(genes.rank, paste(inputName,"_genes.rank",sep = ''),sep = '\t')

dataset.tsne <- data.frame(tsne$Y,lymphocyte.ratio = as.numeric(biospecimen.file), size = as.numeric(xCellScore.brca[65,]))
colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet, lymphocyte.ratio = as.numeric(biospecimen.file), size = as.numeric(xCellScore.brca[65,]))
dataset <- data.frame(example_large_scale$ydata, lymphocyte.ratio = as.numeric(biospecimen.file), size = as.numeric(xCellScore.brca[65,]))
colnames(dataset)[c(1,2)] <- c("component1","component2")

p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,col = lymphocyte.ratio, size = size))
p.tsne <- p.tsne+geom_point()+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))
p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient(low = 'blue', high = 'red')

p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,col = lymphocyte.ratio,size = size))
p.pca <- p.pca+geom_point()+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
p.pca <- p.pca+scale_colour_gradient(low = 'blue', high = 'red')

p.simlr <- ggplot(dataset,aes(x=component1,y=component2,col = lymphocyte.ratio,size = size))
p.simlr <- p.simlr+geom_point()+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))
p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient(low = 'blue', high = 'red')


pdf(file=paste(inputName,"_",featureName,"No_Mutation_Correct.pdf",sep=''),width=28, height=7)
ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
dev.off()




