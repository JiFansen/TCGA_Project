################################################################################################################

### After dN/dS Calculation for each cancre type, we need to bind the mutation table together.

### Then, we need to calculate the score for each possible codon for each significantly mutated gene.

################################################################################################################
cancername <- c("DLBC","GBM","KIRC","LAML","LIHC","OV","PCPG","THYM","UVM","CHOL","ESCA","KICH","KIRP","LGG","MESO","PAAD","TGCT","UCS",
                "CESC","BRCA","HNSC","BLCA","STAD","LUSC","LUAD","COAD","UCEC","THCA","SARC","PRAD","READ","SKCM")
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
####################################################################################################################
strand.plus.missense <- data.initial
uniq.gene <- unique(strand.plus.missense$gene)

happy.results <- data.frame()
genename <- character()
for(i in 1:length(uniq.gene)){
  genename[i] <- as.character(uniq.gene[i])
  test <- strand.plus.missense[which(strand.plus.missense$gene==genename[i]),]
  uniq.codon <- as.character(unique(test$protein.position))
  count <- numeric()
  total.value <- numeric()
  codon.position <- character()
  for(j in 1:length(uniq.codon)){
    codon.position[j] <- uniq.codon[j]
    test.1 <- test[which(test$protein.position==codon.position[j]),]
    count[j] <- dim(test.1)[1]
    total.value[j] <- sum(as.numeric(as.character((test.1$true.value))))
  }
  happy <- data.frame(gene = genename[i],codon.position = codon.position, count = count, total.value = total.value)
  happy.results <- rbind.data.frame(happy.results,happy)
}
happy.results$absolute <- happy.results[,3]/happy.results[,4]
happy.results$`1/absolute` <- 1/happy.results[,5]

score.correction <- data.frame()
for (i in 1:length(uniq.gene)){
  ctnnb1 <- happy.results[which(happy.results$gene==as.character(uniq.gene[i])),]
  ctnnb1$scale <- (ctnnb1$`1/absolute`-min(ctnnb1$`1/absolute`))/(max(ctnnb1$`1/absolute`)-min(ctnnb1$`1/absolute`))
  score.correction <- rbind.data.frame(score.correction,ctnnb1)
}
write.table(score.correction,file="score.correction",sep="\t")