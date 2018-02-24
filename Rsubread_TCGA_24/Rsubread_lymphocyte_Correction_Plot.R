cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","LAML",
                      "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                      "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                      "UCS","UVM","ESCA","GBM","SKCM")
final.biospe.file <- data.frame()
for(i in 1:length(cancername.total)){
  biospecimen.name <- paste(cancername.total[i],".clean.no.zero")
  biospecimen.file <- read.table(file = biospecimen.name, header = T, sep = "\t")
  biospecimen.file <- biospecimen.file[8,]
  biospecimen.file <- t(biospecimen.file)
  colnames(biospecimen.file) <- "lymphocyte_infiltration"
  final.biospe.file <- rbind.data.frame(final.biospe.file,biospecimen.file)
}
final.biospe.file$sample.names <- rownames(final.biospe.file)

data.rsubread24 <- read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt",header=T,sep="\t",rownames=1)
label.biospe <- substr(as.character(rownames(final.biospe.file)),1,16)
label.rsubread <- substr(as.character(colnames(data.rsubread24)),1,16)
rsubread_2_biospe <- match(label.rsubread,label.biospe)
index.rsubread <- which(is.na(rsubread_2_biospe)==FALSE)
data.rsubread.subset <- data.rsubread24[,index.rsubread]
index.biospe <- rsubread_2_biospe[index.rsubread]
data.biospe.subset <- final.biospe.file[index.biospe,]
data.biospe.subset[,1] <- as.numeric(data.biospe.subset[,1])/100
