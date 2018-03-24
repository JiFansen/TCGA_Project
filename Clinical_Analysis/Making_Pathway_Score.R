library('GSVA')
library('preprocessCore')


data.rsubread <- read.table(file = 'GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt',header = T, sep = '\t', row.names = 1)
data.no.small <- data.rsubread[which(rowSums(data.rsubread<1)!=dim(data.rsubread)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName

antiGen <- read.table(file = 'Tcell-Activation_GeneList.txt', header = T, sep = '\t')
antiGen <- as.list(antiGen)

es.max <- gsva(data.no.small, antiGen, mx.diff=FALSE, verbose=FALSE, parallel.sz=1, method = "ssgsea")
rownames(es.max) <- "Tcell-Activation"
write.table(es.max,"Tcell-ActivationScore.txt",sep='\t')

antiGen <- read.table(file = 'Bcell-Activation_GeneList.txt', header = T, sep = '\t')
antiGen <- as.list(antiGen)

es.max <- gsva(data.no.small, antiGen, mx.diff=FALSE, verbose=FALSE, parallel.sz=1, method = "ssgsea")
rownames(es.max) <- "Bcell-Activation"
write.table(es.max,"Bcell-ActivationScore.txt",sep='\t')
