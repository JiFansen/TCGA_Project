
######################################################################################
##################该段代码用来根据片子信息对表达谱数据进行校正######################
######################################################################################

cancername.total <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","HNSC","SKCM",
                "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
                "PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC",
                "UCS","UVM")
for (i in 1:length(cancername.total)){
  setwd("D:/senior_year/R/SVR")
  cancername <- cancername.total[i]
  cancertype <- paste(cancername,"-tumor_symb.exp",sep="")
  non.repeat.genename <- read.table(file="non-repeatgene.txt",header=T)
  non.repeat.genename <- as.matrix(non.repeat.genename)
  BRCA <- read.table(file=cancertype,sep="\t",header=T)
  BRCA <- BRCA[,-1]
  rownames(BRCA) <- non.repeat.genename
  colnames(BRCA) <- substr(colnames(BRCA),1,21)
  
  setwd("D:/senior_year/R/SVR/sample_sldes_yes/")
  sample.slides.yes <- read.table(file = paste(cancername,"-yes",sep = ""),sep = "\t")
  sample.slides.yes.name <- rownames(sample.slides.yes)
  sample.slides.yes.name <- paste(cancername,".",sample.slides.yes.name,sep = "")
  sample.slides.yes.value <- sample.slides.yes[,1]
  
  #setwd("D:/senior_year/R/SVR/sample_slides_no/")
  #sample.slides.no.name <- read.table(file = paste(cancername,"-samples.slides.no",sep = ""),sep = "\t")
  #sample.slides.no.name <- as.character(sample.slides.no.name[,1])
  #sample.slides.no.name <- substr(sample.slides.no.name,1,21)
  #setwd("D:/senior_year/R/SVR/tunedModelY/")
  #sample.slides.no.value <- read.table(file = paste(cancername,"-tunedModelY",sep = ""),sep = "\t")
  #sample.slides.no.value <- 2^as.numeric(sample.slides.no.value[,1])
  
  setwd("D:/senior_year/R/SVR")
  
  data.yes <- BRCA[,match(sample.slides.yes.name,colnames(BRCA))]
  #data.no <- BRCA[,match(sample.slides.no.name,colnames(BRCA))]
  yesfun <- function(x){
    x/sample.slides.yes.value
  }
  #nofun <- function(x){
  #  x/sample.slides.no.value
  #}
  yes.correct <- t(apply(data.yes,1,yesfun))
  #no.correct <- t(apply(data.no, 1,nofun))
  #data.correct <- cbind.data.frame(yes.correct,no.correct)
  setwd("D:/senior_year/R/SVR/correct/")
  write.table(yes.correct,file = paste(cancername,"_correct",sep = ""),sep = "\t")
  setwd("D:/senior_year/R/SVR")
}
