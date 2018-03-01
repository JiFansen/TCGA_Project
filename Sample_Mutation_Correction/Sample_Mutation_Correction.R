







alias <- read.table(file="GeneAlias",header=T)
map <- function(diff.gene){
  trans <- character()
  for (t in 1:length(diff.gene)){
    tmp.gene <- diff.gene[t]
    count = 0
    for(m in 1:dim(alias)[1]){
      a <- as.character(alias[m,3])
      b <- as.character(alias[m,2])
      if(length(which(unlist(strsplit(a,split = "\\|"))==tmp.gene))>0){
        trans[t] <- as.character(alias[m,2])
        count=count+1
      }else if(diff.gene[t]==b){
        trans[t] <- b
        count=count+1
      }
    }
    if(count==0){
      trans[t] <- NA
    }
  }
  return(trans)
}

for(i in 1:length(cancername.total)){
  setwd("D:/senior_year/R/MAF/correct/")
  correction.file <- read.table(file = paste(cancername.total[i],"correct",sep=""), header = T, sep = "\t")
  setwd("D:/senior_year/R/SVR/")
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
  
  expression.profile <- read.table(file = paste(cancername.total[i],"-tumor_symb.exp",sep=""), header = T, sep = "\t")
  gene <- as.character(expression.profile[,1])
  expression.profile <- expression.profile[,-1]
  
  
  diff.gene <- correction.file[which(is.na(match(as.character(correction.file$gene),gene))),]
  diff.gene <- as.character(diff.gene$gene)
  diff.gene
  trans <- map(diff.gene)
  trans
  
  
  genename.final <- as.character(genename.initial)
  genename.final[which(is.na(match(as.character(correction.file$gene),gene)))] <- trans
  correction.file$gene <- genename.final
  strange.index <- which(is.na(as.character(correction.file$gene)))
  if(length(strange.index)>0){
    correction.file <- correction.file[-strange.index,]
  }
  
  if(nchar(cancername.total[i])==3){
    colname <- substr(colnames(expression.profile),5,20)
  }
  if(nchar(cancername.total[i])==4){
    colname <- substr(colnames(expression.profile),6,21)
  }
  if(nchar(cancername.total[i])==2){
    colname <- substr(colnames(expression.profile),4,19)
  }
  colname <- gsub(pattern = "\\.", replacement = "-", x = colname)
  colnames(expression.profile) <- colname
  common.samples <- intersect(as.character(unique(correction.file$sample)),colname)
  expression.profile <- expression.profile[,match(common.samples,colnames(expression.profile))]
  for(j in 1:length(common.samples)){
    inter.test <- correction.file[which(correction.file$sample==common.samples[j]),]
    express.index <- which(is.na(match(gene,as.character(inter.test$gene)))==FALSE)
    inter.test.index <- match(gene,as.character(inter.test$gene))[express.index]
    expression.profile[express.index,j] <- expression.profile[express.index,j]*inter.test[inter.test.index,3]
  }
  setwd("D:/senior_year/R/SVR/Expression_Correct_Data/")
  write.table(expression.profile,file=paste(cancername.total[i],".Correct",sep=""))
  setwd("D:/senior_year/R/SVR/")
}
