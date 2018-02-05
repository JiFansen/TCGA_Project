
####################################################################################
################把所有癌型整理出来的突变表合并成一个大的数据框####################
#####################################################################################

cancername <- c("DLBC","GBM","KIRC","LAML","LIHC","OV","PCPG","THYM","UVM","CHOL","ESCA","KICH","KIRP","LGG","MESO","PAAD","TGCT","UCS",
                "CESC","BRCA","HNSC","BLCA","STAD","LUSC","LUAD","COAD","UCEC","THCA","SARC","PRAD","READ")
data.initial <- read.table(file="/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/results/ACC",header=T,sep="\t")
colname <- colnames(data.initial)
data.initial <- as.matrix(data.initial)
for(i in 1:length(cancername)){
  cancertype <- cancername[i]
  data.between <- read.table(file=cancertype,header=F,sep="\t")
  data.between <- data.between[,-1]
  data.between <- as.matrix(data.between)
  data.initial <- rbind(data.initial,data.between)
}
data.initial <- as.data.frame(data.initial)
genename <- as.character(data.initial$gene)
table <- as.data.frame(table(genename))
table[which(table$Freq==1),]



for(i in 1:length(cancername)){
  cancertype <- cancername[i]
  data.between <- read.table(file=paste(cancertype,".Correct",sep=""),header=T,sep=" ")
  data.between <- as.matrix(data.between)
  data.initial <- cbind(data.initial,data.between)
}
a=1;b=78;
data.new <- data.no.zero[,1:78]
results <- read.table(file="samplelLabels",header=T,sep="\t")
count <- as.numeric(results$Freq)
for(i in 2:32){
  a=b+1;
  b=a+count[i]-1
  test <- data.no.zero[,a:b]
  test <- t(scale(t(test),center = TRUE, scale = TRUE))
  data.new <- cbind(data.new,test)
}

initial <- data.frame(cancertype=rep("ACC",78))
for (i in 2:dim(results)[1]){
  test <- data.frame(cancertype=rep(as.character(results[i,1]),as.numeric(results[i,2])))
  cancertype <- rbind.data.frame(initial, test)
}




cancername <- c("DLBC","GBM","KIRC","LAML","LIHC","OV","PCPG","THYM","UVM","CHOL","ESCA","KICH","KIRP","LGG","MESO","PAAD","TGCT","UCS",
                "CESC","BRCA","HNSC","BLCA","STAD","LUSC","LUAD","COAD","UCEC","THCA","SARC","PRAD","READ")


cancertype <- c(rep("ACC",78),rep("DLBC",37),rep("GBM",148),rep("KIRC",332),rep("LAML",56),rep("LIHC",356),rep("OV",271),rep("PCPG",149),
                rep("THYM",93),rep("UVM",80),rep("CHOL",35),rep("ESCA",160),rep("KICH",64),rep("KIRP",278),rep("LGG",498),rep("MESO",78),
                rep("PAAD",141),rep("TGCT",135),rep("UCS",56),rep("CESC",285),rep("BRCA",969),rep("HNSC",492),rep("BLCA",407),rep("STAD",365),
                rep("LUSC",486),rep("LUAD",502),rep("COAD",396),rep("UCES",522),rep("THCA",450),rep("SARC",235),rep("PRAD",470),rep("READ",132))

final <- character()
for(i in 1:dim(xcellData)[1]){
  test <- as.character(as.matrix(xcellData[i,]))
  final <- union(final,test)
}



