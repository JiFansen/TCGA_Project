######################Read the Data.##########################################
starttime <- Sys.time()
setwd("D:/senior_year/R/MAF/")
filename <- "TCGA.LUAD.varscan.csv"
cancertype <- unlist(strsplit(filename,"\\."))[2]
CHOH <- read.csv(file=filename,header=T)
#header <- colnames(CHOH)
#colnames(CHOH) <- header
useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                 "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                 "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
data <- CHOH[,match(useful.name,colnames(CHOH))]
data <- data[which(data$VARIANT_CLASS=="SNV"),]
data <- data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
#missense.data <- data[which(data$Variant_Classification==c("Missense_Mutation","Silent")),]
data.1 <- data[which(data$Variant_Classification=="Missense_Mutation"),]
data.2 <- data[which(data$Variant_Classification=="Silent"),]
##########################Making Codon Table and Map them all to the DNA.##########################################################
missense.data <- rbind.data.frame(data.1,data.2)
uniq.missence.genes <- unique(missense.data$Hugo_Symbol)
codon <- read.table(file="./codon1.txt",sep="\t",header=F)
codon <- codon[,3:4]
colnames(codon) <- c("Amino Acid","Codon")
codon.table <- matrix(ncol = 2)
for (i in 1:dim(codon)[1]){
  test <- unlist(strsplit(as.character(codon$Codon[i]),split=","))
  matrix <- as.matrix(data.frame(x=rep(codon$`Amino Acid`[i],length(unlist(strsplit(
    as.character(codon$Codon[i]),split=",")))),y=unlist(strsplit(
      as.character(codon$Codon[i]),split=","))))
  codon.table <- rbind(codon.table,matrix)
}
codon.table <- as.data.frame(codon.table[-1,])
colnames(codon.table) <- c("Amino.Acid","Codon")
codon.table <- codon.table[-62,]

for (i in 1:dim(codon.table)[1]){
  char <- unlist(strsplit(as.character(codon.table$Codon[i]),""))
  for (j in 1:length(char)){
    if(char[j]=="U")
    {
      char[j] = "T"
    }
  }
  codon.table[i,3] <- paste(char[1],char[2],char[3],sep="")
}
colnames(codon.table)[3] <- "Coding.strand"

##########For a certain cancer type, extract every useful information for every possible mutation gene.##########
final <- data.frame()
for (i in 1:length(uniq.missence.genes)){
  genename <- as.character(uniq.missence.genes[i])
  test <- missense.data[which(missense.data$Hugo_Symbol==genename),]
  results.2 <- character();results.3 <- character()
  results.4 <- character();results.5 <- character();results.6 <- character();results.7 <- character();results.8 <- character();results.9 <- character()
  for (j in 1:dim(test)[1]){
    results.2[j] <- as.character(test[j,5])
    results.3[j] <- paste(test[j,2],test[j,3],"(",substr(as.character(test[j,6]),5,5),substr(as.character(test[j,6]),7,7),")",sep="")
    results.4[j] <- as.character(test[j,6])
    forward <- unlist(strsplit(as.character(test[j,6]),""))
    reverse <- c()
    for (z in 1:length(forward)){
      if(forward[z]=="A"){
        reverse[z]="T"
      }
      if(forward[z]=="T"){
        reverse[z]="A"
      }
      if(forward[z]=="C"){
        reverse[z]="G"
      }
      if(forward[z]=="G"){
        reverse[z]="C"
      }
    }
    reverse <- paste(reverse[11],reverse[10],reverse[9],reverse[8],reverse[7],reverse[6],reverse[5],reverse[4],reverse[3],reverse[2],reverse[1],sep="")
    results.5[j] <- reverse
    mutation <- unlist(strsplit(as.character(test[j,5]),""))[which(unlist(strsplit(toupper(as.character(test[j,5])),""))==unlist(strsplit(as.character(test[j,5]),"")))]
    position <- which(unlist(strsplit(toupper(as.character(test[j,5])),""))==unlist(strsplit(as.character(test[j,5]),"")))[1]
    codon.mutation <- paste(mutation[1],mutation[3],sep="")
    results.6[j] <- as.character(codon.mutation)
    results.7[j] <- as.character(position)
    results.8[j] <- as.character(test[j,8])
    results.9[j] <- as.character(test[j,7])
  }
  results <- data.frame(gene = genename, codon = results.2, change = results.3, forward = results.4, reverse = results.5, codon.change = results.6, 
                        codon.position = results.7, varaint.class = results.8, protein.position = results.9)
  final <- rbind.data.frame(final,results)
}

######################Judge whether the codon is on the + strand.#######################################
results.10 <- character()
for (i in 1:dim(final)[1]){
  reference.change <- substr(as.character(final[i,3]),1,2)
  codon.change <- as.character(final[i,6])
  if(reference.change==codon.change){
    results.10[i] <- "+"
  }
  else{
    results.10[i] <- "-"
  }
}
final$strand <- results.10

strand.plus <- final[which(final$strand=="+"),]
strand.subtract <- final[which(final$strand=="-"),]
for(i in 1:dim(strand.plus)[1]){
  if(strand.plus[i,7]==1){
    strand.plus[i,11] <- substr(as.character(strand.plus[i,4]),5,9)
  }
  if(strand.plus[i,7]==2){
    strand.plus[i,11] <- substr(as.character(strand.plus[i,4]),4,8)
  }
  if(strand.plus[i,7]==3){
    strand.plus[i,11] <- substr(as.character(strand.plus[i,4]),3,7)
  }
}
colnames(strand.plus)[11] <- "useful.context"

#strand.subtract <- strand.subtract[-40733,]
for(i in 1:dim(strand.subtract)[1]){
  if(strand.subtract[i,7]==1){
    strand.subtract[i,11] <- substr(as.character(strand.subtract[i,5]),5,9)
  }
  if(strand.subtract[i,7]==2){
    strand.subtract[i,11] <- substr(as.character(strand.subtract[i,5]),4,8)
  }
  if(strand.subtract[i,7]==3){
    strand.subtract[i,11] <- substr(as.character(strand.subtract[i,5]),3,7)
  }
}
colnames(strand.subtract)[11] <- "useful.context"

background <- read.table(file="./background_cancer.txt",header=T)
number <- which(colnames(background)==cancertype)
background <- background[,c(1,2,number)]
a <- unlist(strsplit(as.character(background[,2]),split="\\."))
for (i in 1:192){
  b <- paste("(",a[2*i-1],a[2*i],")",sep = "")
  background[i,4] <- b
}
colnames(background)[4] <- "context"

##################Calculate all the possible mutation conditions for a given codon. ##############################################
strand.plus <- rbind.data.frame(strand.plus,strand.subtract)
for(i in 1:dim(strand.plus)[1]){
  string <- substr(as.character(strand.plus[i,11]),2,4)
  amino.acid <- codon.table[which(as.character(codon.table$Coding.strand)==string),1]
  string <- unlist(strsplit(string,""))
  char <- data.frame()
  if(string[1]=="A"){
    char[1,1]<- paste("T",string[2],string[3],sep="")
    char[1,2]<- "AT"
    char[2,1] <- paste("C",string[2],string[3],sep="")
    char[2,2] <- "AC"
    char[3,1] <- paste("G",string[2],string[3],sep="")
    char[3,2] <- "AG"
  }
  if(string[1]=="T"){
    char[1,1] <- paste("A",string[2],string[3],sep="")
    char[1,2] <- "TA"
    char[2,1] <- paste("C",string[2],string[3],sep="")
    char[2,2] <- "TC"
    char[3,1] <- paste("G",string[2],string[3],sep="")
    char[3,2] <- "TG"
  }
  if(string[1]=="C"){
    char[1,1] <- paste("T",string[2],string[3],sep="")
    char[1,2] <- "CT"
    char[2,1] <- paste("A",string[2],string[3],sep="")
    char[2,2] <- "CA"
    char[3,1] <- paste("G",string[2],string[3],sep="")
    char[3,2] <- "CG"
  }
  if(string[1]=="G"){
    char[1,1] <- paste("T",string[2],string[3],sep="")
    char[1,2] <- "GT"
    char[2,1] <- paste("C",string[2],string[3],sep="")
    char[2,2] <- "GC"
    char[3,1] <- paste("A",string[2],string[3],sep="")
    char[3,2] <- "GA"
  }
  if(string[2]=="A"){
    char[4,1] <- paste(string[1],"T",string[3],sep="")
    char[4,2] <- "AT"
    char[5,1] <- paste(string[1],"C",string[3],sep="")
    char[5,2] <- "AC"
    char[6,1] <- paste(string[1],"G",string[3],sep="")
    char[6,2] <- "AG"
  }
  if(string[2]=="T"){
    char[4,1] <- paste(string[1],"A",string[3],sep="")
    char[4,2] <- "TA"
    char[5,1] <- paste(string[1],"C",string[3],sep="")
    char[5,2] <- "TC"
    char[6,1] <- paste(string[1],"G",string[3],sep="")
    char[6,2] <- "TG"
  }
  if(string[2]=="C"){
    char[4,1] <- paste(string[1],"T",string[3],sep="")
    char[4,2] <- "CT"
    char[5,1] <- paste(string[1],"A",string[3],sep="")
    char[5,2] <- "CA"
    char[6,1] <- paste(string[1],"G",string[3],sep="")
    char[6,2] <- "CG"
  }
  if(string[2]=="G"){
    char[4,1] <- paste(string[1],"T",string[3],sep="")
    char[4,2] <- "GT"
    char[5,1] <- paste(string[1],"C",string[3],sep="")
    char[5,2] <- "GC"
    char[6,1] <- paste(string[1],"A",string[3],sep="")
    char[6,2] <- "GA"
  }
  if(string[3]=="A"){
    char[7,1] <- paste(string[1],string[2],"T",sep="")
    char[7,2] <- "AT"
    char[8,1] <- paste(string[1],string[2],"C",sep="")
    char[8,2] <- "AC"
    char[9,1] <- paste(string[1],string[2],"G",sep="")
    char[9,2] <- "AG"
  }
  if(string[3]=="T"){
    char[7,1] <- paste(string[1],string[2],"A",sep="")
    char[7,2] <- "TA"
    char[8,1] <- paste(string[1],string[2],"C",sep="")
    char[8,2] <- "TC"
    char[9,1] <- paste(string[1],string[2],"G",sep="")
    char[9,2] <- "TG"
  }
  if(string[3]=="C"){
    char[7,1] <- paste(string[1],string[2],"T",sep="")
    char[7,2] <- "CT"
    char[8,1] <- paste(string[1],string[2],"A",sep="")
    char[8,2] <- "CA"
    char[9,1] <- paste(string[1],string[2],"G",sep="")
    char[9,2] <- "CG"
  }
  if(string[3]=="G"){
    char[7,1] <- paste(string[1],string[2],"T",sep="")
    char[7,2] <- "GT"
    char[8,1] <- paste(string[1],string[2],"C",sep="")
    char[8,2] <- "GC"
    char[9,1] <- paste(string[1],string[2],"A",sep="")
    char[9,2] <- "GA"
  }
  
  for(j in 1:9){
    if(codon.table[which(as.character(codon.table[,3])==as.character(char[j,1])),1]==amino.acid){
      char[j,3] <- "sys"
    }
    else{
      char[j,3] <- "non"
    }
  }
  for (z in 1:3){
    char[z,4] <- paste("(",unlist(strsplit(as.character(strand.plus[i,11]),""))[1],unlist(strsplit(as.character(strand.plus[i,11]),""))[3],")",sep="")
  }
  for(z in 4:6){
    char[z,4] <- paste("(",unlist(strsplit(as.character(strand.plus[i,11]),""))[2],unlist(strsplit(as.character(strand.plus[i,11]),""))[4],")",sep="")
  }
  for (z in 7:9){
    char[z,4] <- paste("(",unlist(strsplit(as.character(strand.plus[i,11]),""))[3],unlist(strsplit(as.character(strand.plus[i,11]),""))[5],")",sep="")
  }
  for(j in 1:9){
    for(z in 1:dim(background)[1]){
      if(as.character(char[j,2])==as.character(background[z,1])){
        if(as.character(char[j,4])==as.character(background[z,4])){
          char[j,5] <- background[z,3]
        }
      }
    }
  }
  strand.plus[i,12] <- sum(char[which(char[,3]=="non"),5])
  strand.plus[i,13] <- sum(char[which(char[,3]=="sys"),5])
}
colnames(strand.plus)[12] <- "non.background"
colnames(strand.plus)[13] <- "sys.background"
#write.table(strand.plus,file = "sort1_generesults.txt")

final.results <- data.frame()
uniq.plus.genes <- unique(strand.plus$gene)
for (i in 1:length(uniq.plus.genes)){
  gene.name <- as.character(uniq.plus.genes[i])
  test.data <- strand.plus[which(strand.plus$gene==gene.name),]
  non.count <- length(which(test.data$varaint.class=="Missense_Mutation"))
  sys.count <- length(which(test.data$varaint.class=="Silent"))
  #missense.test.data <- test.data[which(test.data$variant.class=="Missense_Mutation"),]
  #silent.test.data <- test.data[which(test.data$variant.class=="Silent"),]
  uniq.protein.position <- unique(test.data$protein.position)
  non.background <- 0
  sys.background <- 0
  for (j in 1:length(uniq.protein.position)){
    protein.position <- uniq.protein.position[j]
    non.background <- non.background+test.data[which(test.data$protein.position==protein.position)[1],12]
    sys.background <- sys.background+test.data[which(test.data$protein.position==protein.position)[1],13]
  }
  final.results[i,1] <- gene.name
  final.results[i,2] <- non.count
  final.results[i,3] <- sys.count
  final.results[i,4] <- non.background
  final.results[i,5] <- sys.background
}
colnames(final.results) <- c("gene","non.count","sys.count","non.background","sys.background")

final.results[,6] <- final.results[,4]/(final.results[,4]+final.results[,5])
final.results[,7] <- final.results[,5]/(final.results[,4]+final.results[,5])
final.results[,2] <- final.results[,2]+1
final.results[,3] <- final.results[,3]+1
final.results[,8] <- (final.results[,2]/final.results[,4])/(final.results[,3]/final.results[,5])
colnames(final.results)[c(6,7,8)] <- c("non.ratio","sys.ratio","dN/ds")
final.results[,9] <- final.results[,2]+final.results[,3]
for(i in 1:dim(final.results)[1]){
  final.results[i,10] <- binom.test((final.results[i,2]-1),(final.results[i,9]-2),final.results[i,6],alternative = "greater")$p.value
}

colnames(final.results)[c(9,10)] <- c("total.count","p.value")
result.1 <- final.results[which(final.results$`dN/ds`>1),]
length(which(result.1$p.value<=0.05))
result <- result.1[which(result.1$p.value<=0.05),]
#write.table(result,file="sort1_significant_results.txt")

value <- numeric()
for (i in 1:dim(strand.plus)[1]){
  mutation <- substr(as.character(strand.plus[i,3]),1,2)
  context <- substr(as.character(strand.plus[i,3]),3,6)
  for (j in 1:192){
    category <- as.character(background[j,1])
    bagrnd.context <- background[j,4]
    if(mutation==category&&context==bagrnd.context)
    {
      value[i]=background[j,3]
    }
  }
}
strand.plus$true.value <- value
genelist <- read.csv(file="./final_binom_sig_gene.csv",header=F)
use <- data.frame()
for(i in 1:dim(genelist)[1]){
  if(length(which(strand.plus$gene==as.character(genelist[i,])))>0){
    test <- strand.plus[which(strand.plus$gene==as.character(genelist[i,])),]
    use <- rbind.data.frame(use,test)
  }
}
use <- use[which(use$varaint.class=="Missense_Mutation"),]
write.table(use,file=cancertype,sep="\t")
#write.table(result,file=paste(cancertype,"result",sep=""),sep="\t")
endtime <- Sys.time()
endtime - starttime