library('preprocessCore')
#setwd("D:/senior_year/R/SVR/")
cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","LAML",
                      "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                      "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                      "UCS","UVM","ESCA","GBM","SKCM")
data.rsubread24 <- read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt",header=T,sep="\t",row.names=1)
xCellScore <- read.table(file="xCellScore.txt",header=T,sep="\t")
xCellScore <- as.matrix(xCellScore)

cancertype <- read.table(file="cancertype",header=F,sep="\t")
cancername <- unique(as.character(cancertype[,2]))
cancertype$sampleBarcode <- substr(as.character(cancertype$V1),1,16)
data.no.small <- data.rsubread24[which(rowSums(data.rsubread24<1)!=dim(data.rsubread24)[2]),]
data.no.small <- log2(data.no.small+1)
data.no.small <- as.matrix(data.no.small)
sampleName <- colnames(data.no.small)
genename <- rownames(data.no.small)
data.no.small <- normalize.quantiles(data.no.small,copy=TRUE)
rownames(data.no.small) <- genename
colnames(data.no.small) <- sampleName
data.final <- data.no.small
data.subset <- matrix(1,nrow=dim(data.no.small)[1])
xCellScore.subset <- matrix(1, nrow = dim(xCellScore)[1])
CancerLabel <- character()
for(i in match(unique(as.character(cancertype$V2)),cancername.total)){
  
  setwd("/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/correct")
  correction.file <- read.table(file = paste(cancername.total[i],"correct",sep=""), header = T, sep = "\t")
  setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24")
  print(i)
  
  genename.initial <- as.character(correction.file$gene)
  genename.initial[which(genename.initial=="7-Mar")] <- "LDOC1"
  genename.initial[which(genename.initial=="2-Mar")] <- "PEG10"
  genename.initial[which(genename.initial=="8-Mar")] <- "RTL8C"
  genename.initial[which(genename.initial=="6-Mar")] <- "RTL6"
  genename.initial[which(genename.initial=="10-Mar")] <- "MARCH10"
  genename.initial[which(genename.initial=="1-Mar")] <- "Mar1"
  genename.initial[which(genename.initial=="4-Mar")] <- "RTL4"
  genename.initial[which(genename.initial=="11-Mar")] <- "MARCH11"
  genename.initial[which(genename.initial=="8-Sep")] <- "SEPT8"
  genename.initial[which(genename.initial=="3-Sep")] <- "SEPT3"
  genename.initial[which(genename.initial=="6-Sep")] <- "SEPT6"
  genename.initial[which(genename.initial=="2-Sep")] <- "SEPT2"
  genename.initial[which(genename.initial=="7-Sep")] <- "SEPT7"
  genename.initial[which(genename.initial=="4-Sep")] <- "SEPT4"
  genename.initial[which(genename.initial=="5-Sep")] <- "SEPT5"
  genename.initial[which(genename.initial=="1-Sep")] <- "SEPT1"
  genename.initial[which(genename.initial=="10-Sep")] <- "SEPT10"
  genename.initial[which(genename.initial=="11-Sep")] <- "SEPT11"
  genename.initial[which(genename.initial=="12-Sep")] <- "SEPT12"
  genename.initial[which(genename.initial=="15-Sep")] <- "SEP15"
  genename.initial[which(genename.initial=="14-Sep")] <- "SEPT14"
  genename.initial[which(genename.initial=="1-Dec")] <- "DEC1"
  
  correction.file$gene <- genename.initial
  common.samples <- intersect(as.character(unique(correction.file$sample)),as.character(cancertype$sampleBarcode))
  
  expression.profile <- data.no.small
  score.between <- xCellScore
  index <- match(common.samples,as.character(cancertype$sampleBarcode))
  score.between <- xCellScore[,index]
  expression.profile <- expression.profile[,index]
  for(j in 1:length(common.samples)){
    inter.test <- correction.file[which(correction.file$sample==common.samples[j]),]
    express.index <- which(is.na(match(genename,as.character(inter.test$gene)))==FALSE)
    inter.test.index <- match(genename,as.character(inter.test$gene))[express.index]
    expression.profile[express.index,j] <- expression.profile[express.index,j]*inter.test[inter.test.index,3]
  }
  data.final[,index] <- expression.profile
  data.subset <- cbind(data.subset,expression.profile)
  CancerLabel <- c(CancerLabel,rep(cancername.total[i],length(common.samples)))
  xCellScore.subset <- cbind(xCellScore.subset,score.between)
}
data.subset <- data.subset[,-1]
xCellScore.subset <- xCellScore.subset[,-1]
#write.table(data.final,"data.final.correction",sep = '\t')
write.table(data.subset,"data.subset",sep = '\t')
write.table(CancerLabel,"CancerLabel",sep = '\t')
write.table(xCellScore.subset,"xCellScore.subset",sep = '\t')
