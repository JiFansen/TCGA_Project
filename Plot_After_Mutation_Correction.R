#######################################################################################

### After mutation correction, we need to bind the correction data of Each cancer type.

### Also, after this step, we need to remove genes that FPKM values are less than 1 in all
### samples to make our data cleaner.

### Then, we need to log2(data+1).

### Then, we should quantile normalization of all cancer types data.

### Next, we should scale the data according to each cancer type.

### After these steps, we can do the following clustering analysis.

### It should be noted that we should record the gene name in each step, maybe some of the 
### gene names are duplicate.

#######################################################################################
cancerName <- c("DLBC","GBM","KIRC","LAML","LIHC","OV","PCPG","THYM","UVM","CHOL","ESCA","KICH","KIRP","LGG","MESO","PAAD","TGCT","UCS",
                "CESC","BRCA","HNSC","BLCA","STAD","LUSC","LUAD","COAD","UCEC","THCA","SARC","PRAD","READ")
data.initial <- read.table(file="ACC.Correct",header=T,sep=" ")
for(i in 1:length(cancerName)){
  cancertype <- cancerName[i]
  data.between <- read.table(file=paste(cancertype,".Correct",sep=""),header=T,sep=" ")
  data.initial <- cbind.data.frame(data.initial,data.between)
}

write.table(data.initial,"dataAfterMutationCorrection",sep="\t")
###########################################################################################
library('preprocessCore')
library(SIMLR)
library(igraph)
library(ggplot2)
data.initial <- read.table(file = "dataAfterMutationCorrection",header = T, sep="\t")
numberOfEachSample <- c(rep("ACC",78),rep("DLBC",37),rep("GBM",148),rep("KIRC",332),rep("LAML",56),rep("LIHC",356),rep("OV",271),rep("PCPG",149),
                        rep("THYM",93),rep("UVM",80),rep("CHOL",35),rep("ESCA",160),rep("KICH",64),rep("KIRP",278),rep("LGG",498),rep("MESO",78),
                        rep("PAAD",141),rep("TGCT",135),rep("UCS",56),rep("CESC",285),rep("BRCA",969),rep("HNSC",492),rep("BLCA",407),rep("STAD",365),
                        rep("LUSC",486),rep("LUAD",502),rep("COAD",396),rep("UCES",522),rep("THCA",450),rep("SARC",235),rep("PRAD",470),rep("READ",132))
nonrepeat.gene <- read.table(file="non-repeatgene.txt",header=T,sep="\t")
nonrepeat.gene <-as.character(nonrepeat.gene[,1])
rownames(data.initial) <- nonrepeat.gene

repeat.gene <- read.table(file="repeatgene.txt",sep="\t",header=T)
repeat.gene <- as.character(repeat.gene[,1])

##################################################################################################
no.small.index <- which(rowSums(data.initial<1)!=dim(data.initial)[2])
data.no.small <- data.initial[no.small.index,]
repeat.gene <- repeat.gene[no.small.index]
##################################################################################################
### In this part, we can replace xCellGeneList to ang geneList.
##################################################################################################
xcellGeneList <- read.table(file = "xCellGeneList",header=T,sep = "\t")
xcellGeneList <- as.character(xcellGeneList[,1])
xcellGeneList <- intersect(xcellGeneList,repeat.gene)
xcell.index <- which(is.na(match(repeat.gene,xcellGeneList))==FALSE)

data.xcell <- data.no.small[xcell.index,]
repeat.gene.xcell <- repeat.gene[xcell.index,]
###################################################################################################
data.xcell <- log2(data.xcell+1)
data.xcell <- as.matrix(data.xcell)
sampleName <- colnames(data.xcell)
data.xcell <- normalize.quantiles(data.xcell,copy = TRUE)
###################################################################################################
a=1;b=78;
data.new <- data.xcell[,1:78]
results <- read.table(file="samplelLabels",header=T,sep="\t")
count <- as.numeric(results$Freq)
for(i in 2:32){
  a=b+1;
  b=a+count[i]-1
  test <- data.xcell[,a:b]
  test <- t(scale(t(test),center = TRUE, scale = TRUE))
  data.new <- cbind(data.new,test)
}
#####################################################################################################
example_large_scale = SIMLR_Large_Scale(X = data.new, c = 8, kk = 10)
dataset <- data.frame(example_large_scale$ydata,cancertype=numberOfEachSample)
colnames(dataset) <- c("component1","component2","cancertype")
png(file="SIMLR.png",width=300, height=100, units='mm',res=1500)
ggplot(dataset,aes(x=component1,y=component2,col=cancertype))+geom_point()
dev.off()
######################################################################################################