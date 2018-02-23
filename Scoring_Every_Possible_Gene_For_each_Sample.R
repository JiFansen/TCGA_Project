#####################################################################################################

#### In this step, we need to calculate the gene score for each sample in a certain cancertype.

#####################################################################################################
setwd("D:/senior_year/R/MAF")
score.correction <- read.table(file="score.correction",header=T,sep="\t")
cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","LAML",
                      "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                      "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                      "UCS","UVM","ESCA","GBM")
for (x in 1:length(cancername.total)){
  filename <- paste("TCGA.",cancername.total[x],".varscan.csv",sep="")
  cancertype <- unlist(strsplit(filename,"\\."))[2]
  setwd("D:/senior_year/R/MAF/MAFile/")
  maf.data <- read.csv(file=filename,header=T)
  setwd("D:/senior_year/R/MAF")
  useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                   "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                   "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
  maf.data <- maf.data[,match(useful.name,colnames(maf.data))]
  maf.data <- maf.data[which(maf.data$VARIANT_CLASS=="SNV"),]
  maf.data <- maf.data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
  maf.data <- maf.data[which(maf.data$Variant_Classification=="Missense_Mutation"),]
  maf.data$Sample_Barcode <- substr(as.character(maf.data$Tumor_Sample_Barcode),1,16)
  maf.data <- maf.data[,-12]
  
  uniq.gene <- unique(as.character(as.matrix(score.correction$gene)))
  uniq.sample <- unique(as.character(as.matrix(maf.data$Sample_Barcode)))
  inter <- numeric()
  for (i in 1:length(uniq.sample)){
    test <- maf.data[which(maf.data$Sample_Barcode==uniq.sample[i]),]
    genename <- as.character(as.matrix(test$Hugo_Symbol))
    inter[i] <- length(intersect(uniq.gene,genename))
  }
  
  condition <- (inter>0)
  result.score <- data.frame()
  for (i in (1:length(uniq.sample))[condition]){
    test <- maf.data[which(maf.data$Sample_Barcode==uniq.sample[i]),]
    genename <- as.character(as.matrix(test$Hugo_Symbol))
    test.inter <- subset(test, Hugo_Symbol%in%intersect(uniq.gene,genename),select = c("Hugo_Symbol","Codons","Protein_position","Sample_Barcode"))
    subset <- numeric()
    for(j in 1:dim(test.inter)[1]){
      subset[j] <- score.correction[which(as.character(score.correction$gene)==test.inter[j,1]&as.character(score.correction$codon.position)==test.inter[j,3]),7]
    }
    test.inter$score <- subset
    result.score <- rbind.data.frame(result.score,test.inter)
  }
  
  final.score <- data.frame()
  for(i in 1:length(unique(result.score$Sample_Barcode))){
    subdata.1 <- result.score[which(as.character(result.score$Sample_Barcode)==unique(result.score$Sample_Barcode)[i]),]
    gene.uniq <- as.character(unique(subdata.1$Hugo_Symbol))
    total.score <- numeric()
    for(j in 1:length(gene.uniq)){
      subdata.2 <- subdata.1[which(as.character(subdata.1$Hugo_Symbol)==gene.uniq[j]),]
      score = 1
      for(z in 1:dim(subdata.2)[1]){
        score = score*subdata.2[z,5]
      }
      total.score[j] <- score
    }
    total.result <- data.frame(gene = gene.uniq, sample = uniq.sample[i], total.score = total.score)
    final.score <- rbind.data.frame(final.score,total.result)
  }
  
  setwd("D:/senior_year/R/MAF/FrameshiftResults/")
  frameshift <- read.table(file=paste(cancertype,"result",sep=""),header=T,sep=" ")
  setwd("D:/senior_year/R/MAF")
  frameshift <- frameshift[,c(2,1)]
  frameshift$total.score <- 0
  colnames(frameshift) <- colnames(final.score)
  mutation.frameshift <- rbind.data.frame(final.score,frameshift)
  setwd("D:/senior_year/R/MAF/mutation.frameshift/")
  write.table(mutation.frameshift,file = paste(cancertype,"mutation.frameshift",sep=""),sep="\t")
  setwd("D:/senior_year/R/MAF")
  sample.uniq <- unique(as.character(mutation.frameshift$sample))
  
  results.2 <- data.frame()
  for (i in 1:length(sample.uniq)){
    subdata.1 <- mutation.frameshift[which(mutation.frameshift$sample==sample.uniq[i]),]
    gene.uniq <- unique(as.character(subdata.1$gene))
    score.final <- as.numeric()
    for (j in 1:length(gene.uniq)){
      subdata.2 <- subdata.1[which(subdata.1$gene==gene.uniq[j]),]
      score.total <- 1
      for (z in 1:dim(subdata.2)[1]){
        score.total <- score.total*subdata.2[z,3]
      }
      score.final[j] <- score.total
    }
    results.1 <- data.frame(gene = gene.uniq, sample = sample.uniq[i],score = score.final)
    results.2 <- rbind.data.frame(results.2,results.1)
  }
  setwd("D:/senior_year/R/MAF/correct/")
  write.table(results.2,file = paste(cancertype,"correct",sep=""),sep="\t")
  setwd("D:/senior_year/R/MAF")
}
