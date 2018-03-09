
setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/")
cancertype <- read.table(file = 'CancerLabel',header = T, sep = '\t')
cancertype <- as.character(cancertype[,1])

cancer.total <- as.character(unique(cancertype))

antiGen.correction.file <- matrix(1, nrow = 77)
antiGen.no.correction.file <- matrix(1, nrow = 77)
xCell.correction.file <- matrix(1, nrow = 2367)
xCell.no.correction.file <- matrix(1, nrow = 2367)

setwd("/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24/Gene_Rank")

for(n in 1:length(cancer.total)){
  antiGen.correction.name <- paste(cancer.total[n], "AntigenProcessing_KEGG.txt_genes.correct.rank", sep = '')
  antiGen.no.correction.name <- paste(cancer.total[n], "AntigenProcessing_KEGG.txt_genes.rank", sep = '')
  xCell.correction.name <- paste(cancer.total[n], "xCellGeneList_genes.correct.rank", sep = '')
  xCell.no.correction.name <- paste(cancer.total[n], "xCellGeneList_genes.rank", sep = '')
  
  data.between <- read.table(file = antiGen.correction.name, header = T, sep = '\t')
  data.between <- as.matrix(data.between[,1])
  colnames(data.between)[1] <- cancer.total[n]
  antiGen.correction.file <- cbind(antiGen.correction.file, data.between)
  
  data.between <- read.table(file = antiGen.no.correction.name, header = T, sep = '\t')
  data.between <- as.matrix(data.between[,1])
  colnames(data.between)[1] <- cancer.total[n]
  antiGen.no.correction.file <- cbind(antiGen.no.correction.file, data.between)
  
  data.between <- read.table(file = xCell.correction.name, header = T, sep = '\t')
  data.between <- as.matrix(data.between[,1])
  colnames(data.between)[1] <- cancer.total[n]
  xCell.correction.file <- cbind(xCell.correction.file, data.between)
  
  data.between <- read.table(file = xCell.no.correction.name, header = T, sep = '\t')
  data.between <- as.matrix(data.between[,1])
  colnames(data.between)[1] <- cancer.total[n]
  xCell.no.correction.file <- cbind(xCell.no.correction.file, data.between)
}

antiGen.correction.file <- antiGen.correction.file[,-1]
antiGen.no.correction.file <- antiGen.no.correction.file[,-1]
xCell.correction.file <- xCell.correction.file[,-1]
xCell.no.correction.file <- xCell.no.correction.file[,-1]


write.table(antiGen.correction.file,"antiGen.correction.file",sep = '\t')
write.table(antiGen.no.correction.file,"antiGen.no.correction.file",sep = '\t')
write.table(xCell.correction.file,"xCell.correction.file",sep = '\t')
write.table(xCell.no.correction.file,"xCell.no.correction.file",sep = '\t')


Results_antiGen.correction <- list()
Results_antiGen.no.correction <- list()
Results_xCell.correction <- list()
Results_xCell.no.correction <- list()
for(m in seq(10, 70, by = 10)){
  genename <- c()
  a <- antiGen.correction.file[1:m,]
  for(i in 1:dim(a)[2]){
    genename <- c(genename,as.character(a[,i]))
  }
  between <- as.data.frame(table(genename))
  between <- between[order(-between$Freq),]
  colnames(between)[2] <- paste("top",m,sep = '')
  Results_antiGen.correction[[m/10]] <- between
}

for(m in seq(10, 70, by = 10)){
  genename <- c()
  a <- antiGen.no.correction.file[1:m,]
  for(i in 1:dim(a)[2]){
    genename <- c(genename,as.character(a[,i]))
  }
  between <- as.data.frame(table(genename))
  between <- between[order(-between$Freq),]
  colnames(between)[2] <- paste("top",m,sep = '')
  Results_antiGen.no.correction[[m/10]] <- between
}

y=1
for(m in c(10,20,50,100,500,1000,1500,2000)){
  genename <- c()
  a <- xCell.correction.file[1:m,]
  for(i in 1:dim(a)[2]){
    genename <- c(genename,as.character(a[,i]))
  }
  between <- as.data.frame(table(genename))
  between <- between[order(-between$Freq),]
  colnames(between)[2] <- paste("top",m,sep = '')
  Results_xCell.correction[[y]] <- between
  y <- y+1
}

y=1
for(m in c(10,20,50,100,500,1000,1500,2000)){
  genename <- c()
  a <- xCell.no.correction.file[1:m,]
  for(i in 1:dim(a)[2]){
    genename <- c(genename,as.character(a[,i]))
  }
  between <- as.data.frame(table(genename))
  between <- between[order(-between$Freq),]
  colnames(between)[2] <- paste("top",m,sep = '')
  Results_xCell.no.correction[[y]] <- between
  y <- y+1
}

xCell.correction.file
antiGen.correction.file

SASH3.cancer <- c();WAS.cancer <- c();NCKAP1L.cancer <- c();CD37.cancer <- c();CD48.cancer <- c()
CD53.cancer <- c(); CYTH4.cancer <- c();AIF1.cancer <- c();DOK2.cancer <- c();IL10RA.cancer <- c()
for(i in 1:dim(xCell.correction.file)[2]){
  a <- xCell.correction.file[1:50,]
  gene.list <- as.character(a[,i])
  if(length(which(gene.list=="SASH3"))==1){
    SASH3.cancer <- c(SASH3.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="WAS"))==1){
    WAS.cancer <- c(WAS.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="NCKAP1L"))==1){
    NCKAP1L.cancer <- c(NCKAP1L.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="CD37"))==1){
    CD37.cancer <- c(CD37.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="CD48"))==1){
    CD48.cancer <- c(CD48.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="CD53"))==1){
    CD53.cancer <- c(CD53.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="CYTH4"))==1){
    CYTH4.cancer <- c(CYTH4.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="AIF1"))==1){
    AIF1.cancer <- c(AIF1.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="DOK2"))==1){
    DOK2.cancer <- c(DOK2.cancer, colnames(a)[i])
  }
  if(length(which(gene.list=="IL10RA"))==1){
    IL10RA.cancer <- c(IL10RA.cancer, colnames(a)[i])
  }
}

HLA.DRA.cancer <- c();HLA.DPA1.cancer <- c();HLA.DRB1.cancer <- c();HLA.DPB1.cancer <- c();CD74.cancer <- c();
HLA.DMA.cancer <- c();HLA.DMB.cancer <- c();HLA.DQA1.cancer <- c(); CIITA.cancer <- c(); HLA.DRB5.cancer <- c()

for(i in 1:dim(antiGen.correction.file)[2]){
  a <- antiGen.correction.file[1:10,]
  gene.list <- as.character(a[,i])
  if(length(which(gene.list=='HLA-DRA'))==1){
    HLA.DRA.cancer <- c(HLA.DRA.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DPA1'))==1){
    HLA.DPA1.cancer <- c(HLA.DPA1.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DRB1'))==1){
    HLA.DRB1.cancer <- c(HLA.DRB1.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DPB1'))==1){
    HLA.DPB1.cancer <- c(HLA.DPB1.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='CD74'))==1){
    CD74.cancer <- c(CD74.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DMA'))==1){
    HLA.DMA.cancer <- c(HLA.DMA.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DMB'))==1){
    HLA.DMB.cancer <- c(HLA.DMB.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DQA1'))==1){
    HLA.DQA1.cancer <- c(HLA.DQA1.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='CIITA'))==1){
    CIITA.cancer <- c(CIITA.cancer,colnames(a)[i])
  }
  if(length(which(gene.list=='HLA-DRB5'))==1){
    HLA.DRB5.cancer <- c(HLA.DRB5.cancer,colnames(a)[i])
  }
}





