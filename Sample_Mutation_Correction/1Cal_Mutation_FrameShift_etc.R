setwd("D:/senior_year/R/MAF")
cancername <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
                "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                "UCS","UVM","SKCM")
for(z in 1:32){
  cancertype <- cancername[z]
  filename <- paste("TCGA.",cancertype,".varscan.csv",sep="")
CHOH <- read.csv(file=filename,header=T)
useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                 "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                 "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
data <- CHOH[,match(useful.name,colnames(CHOH))]
data <- data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
data$tumor.sample <- substr(as.character(data$Tumor_Sample_Barcode),1,16)
data.1 <- data[which(data$Variant_Classification=="Frame_Shift_Del"),]
data.2 <- data[which(data$Variant_Classification=="Splice_Site"),]
data.3 <- data[which(data$Variant_Classification=="Nonsense_Mutation"),]
data.4 <- data[which(data$Variant_Classification=="Frame_Shift_Ins"),]
data.5 <- data[which(data$Variant_Classification=="Splice_Region"),]
#data.6 <- data[which(data$Variant_Classification=="In_Frame_Ins"),]
#data.7 <- data[which(data$Variant_Classification=="In_Frame_Del"),]
data.8 <- data[which(data$Variant_Classification=="Nonstop_Mutation"),]
data.9 <- data[which(data$Variant_Classification=="Translation_Start_Site"),]
data.total <- rbind.data.frame(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
rm(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
result <- data.frame()
for(i in 1:length(unique(data.total$tumor.sample))){
  sample.name <- unique(data.total$tumor.sample)[i]
  test <- data.total[which(data.total$tumor.sample==sample.name),]
  count.1 <- numeric();count.2 <- numeric();count.3 <- numeric();count.4 <- numeric();count.5 <- numeric();
  count.6 <- numeric();count.7 <- numeric();count.8 <- numeric();count.9 <- numeric();gene.name <- character()
  for(j in 1:length(unique(test$Hugo_Symbol))){
    gene.name[j] <- as.character(unique(test$Hugo_Symbol)[j])
    test.1 <- test[which(test$Hugo_Symbol==gene.name[j]),]
    count.1[j] <- length(which(test.1$Variant_Classification=="Frame_Shift_Del"))
    count.2[j] <- length(which(test.1$Variant_Classification=="Splice_Site"))
    count.3[j] <- length(which(test.1$Variant_Classification=="Nonsense_Mutation"))
    count.4[j] <- length(which(test.1$Variant_Classification=="Frame_Shift_Ins"))
    count.5[j] <- length(which(test.1$Variant_Classification=="Splice_Region"))
    #count.6[j] <- length(which(test.1$Variant_Classification=="In_Frame_Ins"))
    #count.7[j] <- length(which(test.1$Variant_Classification=="In_Frame_Del"))
    count.8[j] <- length(which(test.1$Variant_Classification=="Nonstop_Mutation"))
    count.9[j] <- length(which(test.1$Variant_Classification=="Translation_Start_Site"))
  }
  result.1 <- data.frame(SampleName = sample.name, GeneName = gene.name, Frame.Shift.Del = count.1, Splice.Site = count.2, Nonsense.Mutation = count.3, Frame.Shift.Ins = count.4, 
                         Splice.Region = count.5,  Nonstop.Mutation = count.8, Translation.Start.Site = count.9)
  result <- rbind.data.frame(result, result.1)
}
count.sum <- numeric()
for(i in 1:dim(result)[1]){
  count.sum[i] <- sum(result[i,c(3:9)])
}
result$sum <- count.sum
final <- paste(cancertype,"result",sep="")
setwd("D:/senior_year/R/MAF/FrameshiftResults/")
write.table(result,file=final)
setwd("D:/senior_year/R/MAF")
}
z=33
cancertype <- cancername[z]
CHOH <- read.csv(file="TCGA.SKCM.varscan.txt",header=T,sep = "\t")
useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                 "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                 "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
data <- CHOH[,match(useful.name,colnames(CHOH))]
data <- data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
data$tumor.sample <- substr(as.character(data$Tumor_Sample_Barcode),1,16)
data.1 <- data[which(data$Variant_Classification=="Frame_Shift_Del"),]
data.2 <- data[which(data$Variant_Classification=="Splice_Site"),]
data.3 <- data[which(data$Variant_Classification=="Nonsense_Mutation"),]
data.4 <- data[which(data$Variant_Classification=="Frame_Shift_Ins"),]
data.5 <- data[which(data$Variant_Classification=="Splice_Region"),]
#data.6 <- data[which(data$Variant_Classification=="In_Frame_Ins"),]
#data.7 <- data[which(data$Variant_Classification=="In_Frame_Del"),]
data.8 <- data[which(data$Variant_Classification=="Nonstop_Mutation"),]
data.9 <- data[which(data$Variant_Classification=="Translation_Start_Site"),]
data.total <- rbind.data.frame(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
rm(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
result <- data.frame()
for(i in 1:length(unique(data.total$tumor.sample))){
  sample.name <- unique(data.total$tumor.sample)[i]
  test <- data.total[which(data.total$tumor.sample==sample.name),]
  count.1 <- numeric();count.2 <- numeric();count.3 <- numeric();count.4 <- numeric();count.5 <- numeric();
  count.6 <- numeric();count.7 <- numeric();count.8 <- numeric();count.9 <- numeric();gene.name <- character()
  for(j in 1:length(unique(test$Hugo_Symbol))){
    gene.name[j] <- as.character(unique(test$Hugo_Symbol)[j])
    test.1 <- test[which(test$Hugo_Symbol==gene.name[j]),]
    count.1[j] <- length(which(test.1$Variant_Classification=="Frame_Shift_Del"))
    count.2[j] <- length(which(test.1$Variant_Classification=="Splice_Site"))
    count.3[j] <- length(which(test.1$Variant_Classification=="Nonsense_Mutation"))
    count.4[j] <- length(which(test.1$Variant_Classification=="Frame_Shift_Ins"))
    count.5[j] <- length(which(test.1$Variant_Classification=="Splice_Region"))
    #count.6[j] <- length(which(test.1$Variant_Classification=="In_Frame_Ins"))
    #count.7[j] <- length(which(test.1$Variant_Classification=="In_Frame_Del"))
    count.8[j] <- length(which(test.1$Variant_Classification=="Nonstop_Mutation"))
    count.9[j] <- length(which(test.1$Variant_Classification=="Translation_Start_Site"))
  }
  result.1 <- data.frame(SampleName = sample.name, GeneName = gene.name, Frame.Shift.Del = count.1, Splice.Site = count.2, Nonsense.Mutation = count.3, Frame.Shift.Ins = count.4, 
                         Splice.Region = count.5,  Nonstop.Mutation = count.8, Translation.Start.Site = count.9)
  result <- rbind.data.frame(result, result.1)
}
count.sum <- numeric()
for(i in 1:dim(result)[1]){
  count.sum[i] <- sum(result[i,c(3:9)])
}
result$sum <- count.sum
final <- paste(cancertype,"result",sep="")
setwd("D:/senior_year/R/MAF/FrameshiftResults/")
write.table(result,file=final)