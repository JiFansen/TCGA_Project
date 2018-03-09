library('preprocessCore')
library("SIMLR")
library("igraph")
library("ggplot2")
library("Rtsne")
library("ggpubr")

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

cancertype <- read.table(file = 'cancertype',header = F, sep = '\t')
cancertype <- as.character(cancertype[,2])

xCellScore <- read.table(file="xCellScore.txt",header=T,sep="\t")

xCellGene <- read.table(file=featureName,header=T,sep="\t")
xCellGene <- as.character(xCellGene[,1])
xcellGeneList <- intersect(xCellGene,genename)
xcell.index <- which(is.na(match(genename,xcellGeneList))==FALSE)
data.xcell <- data.no.small[xcell.index,]

cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","LAML",
                      "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                      "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                      "UCS","UVM","ESCA","GBM","SKCM")
final.biospe.file <- data.frame()
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Slides_Data")
for(i in 1:length(cancername.total)){
  biospecimen.name <- paste(cancername.total[i],".clean.no.zero")
  biospecimen.file <- read.table(file = biospecimen.name, header = T, sep = "\t")
  biospecimen.file <- biospecimen.file[8,]
  biospecimen.file <- t(biospecimen.file)
  colnames(biospecimen.file) <- "lymphocyte_infiltration"
  final.biospe.file <- rbind.data.frame(final.biospe.file,biospecimen.file)
}
final.biospe.file$sample.names <- rownames(final.biospe.file)
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")

expression.sample.name <- substr(colnames(data.xcell),1,16)
biospeci.sample.name <- substr(rownames(final.biospe.file),1,16)

match.vector <- match(expression.sample.name,biospeci.sample.name)
expression.index <- which(is.na(match.vector)==FALSE)
biospeci.index <- match.vector[expression.index]

data.xcell <- data.xcell[,expression.index]
cancertype <- cancertype[expression.index]
final.biospe.file <- final.biospe.file[biospeci.index,]
xCellScore.subset <- xCellScore[,expression.index]

data.xcell <- t(scale(t(data.xcell),center = TRUE, scale = TRUE))
cancername <- unique(cancertype)
#for(i in 1:length(cancername)){
#  cancer <- cancername[i]
#  samples.Name <- which(cancertype==cancer)
#  data.between <- data.xcell[,samples.Name]
#  data.between <- t(scale(t(data.between),center = TRUE, scale = TRUE))
#  data.xcell[,samples.Name] <- data.between
#}

cancershape <- c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','%','h','@','&','#','$','^')

tsne <- Rtsne(t(data.xcell), dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
data.pca <- t(data.xcell)
data.pr <- prcomp(data.pca)
example_large_scale = SIMLR(X = data.xcell, c = 8)

dataset.tsne <- data.frame(tsne$Y,lymphocyte.ratio = as.numeric(final.biospe.file[,1]), size = as.numeric(xCellScore.subset[65,]), cancertype = cancertype)
colnames(dataset.tsne)[c(1,2)] <- c("component1","component2")
dataSet <- data.pr$x[,1:2]
dataSet <- data.frame(dataSet, lymphocyte.ratio = as.numeric(final.biospe.file[,1]), size = as.numeric(xCellScore.subset[65,]), cancertype = cancertype)
dataset <- data.frame(example_large_scale$ydata, lymphocyte.ratio = as.numeric(final.biospe.file[,1]), size = as.numeric(xCellScore.subset[65,]), cancertype = cancertype)
colnames(dataset)[c(1,2)] <- c("component1","component2")


p.tsne <- ggplot(dataset.tsne,aes(x=component1,y=component2,col = lymphocyte.ratio, size = size, shape = cancertype))
p.tsne <- p.tsne+geom_point()+labs(title="tsne")+theme(plot.title=element_text(hjust=0.5))
p.tsne <- p.tsne+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient(low = 'blue', high = 'red')
p.tsne <- p.tsne+scale_shape_manual(values=cancershape[1:length(cancername)])

p.pca <- ggplot(dataSet,aes(x=PC1,y=PC2,col = lymphocyte.ratio,size = size, shape = cancertype))
p.pca <- p.pca+geom_point()+labs(title="pca")+theme(plot.title=element_text(hjust=0.5))+scale_size(range=c(0,2))
p.pca <- p.pca+scale_colour_gradient(low = 'blue', high = 'red')
p.pca <- p.pca+scale_shape_manual(values=cancershape[1:length(cancername)])

p.simlr <- ggplot(dataset,aes(x=component1,y=component2,col = lymphocyte.ratio,size = size, shape = cancertype))
p.simlr <- p.simlr+geom_point()+labs(title="simlr")+theme(plot.title=element_text(hjust=0.5))
p.simlr <- p.simlr+scale_size(range=c(0,2))+theme(legend.position="none")+scale_colour_gradient(low = 'blue', high = 'red')
p.simlr <- p.simlr+scale_shape_manual(values=cancershape[1:length(cancername)])


pdf(file=paste(featureName,"All_Cancertype_No_Mutation_Correct.pdf",sep=''),width=28, height=7)
ggarrange(p.tsne, p.pca, p.simlr,ncol=3,nrow=1,labels=c("A","B","C"))
dev.off()










