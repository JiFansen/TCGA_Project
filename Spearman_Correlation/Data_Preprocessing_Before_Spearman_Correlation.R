############################################################################################

### At this step, we need to read in the expression matrix for each cancer type and the 

### lymphocyte slides data. After normalize the expression data, we need to match the samples

### names of these two types of data. This part will tell us how many samples are duplicated?

##############################################################################################
cancername <- "ESCA"
cancertype <- paste(cancername,"-tumor_symb.exp",sep="")
non.repeat.genename <- read.table(file="non-repeatgene.txt",header=T)
non.repeat.genename <- as.matrix(non.repeat.genename)
BRCA <- read.table(file=cancertype,sep="\t",header=T)
BRCA <- BRCA[,-1]
rownames(BRCA) <- non.repeat.genename
BRCA <- BRCA[which(rowSums(BRCA>=1)!=0),]
gene <- rownames(BRCA)
sample <- colnames(BRCA)
BRCA <- as.matrix(BRCA)
library(preprocessCore)
BRCA <- log2(BRCA+1)
BRCA <- normalize.quantiles(BRCA,copy=TRUE)
BRCA <- t(scale(t(BRCA), center=T, scale = T))
colnames(BRCA) <- sample
rownames(BRCA) <- gene
lymphocyte.name <- paste(cancername," .clean.no.zero",sep="")
slides.BRCA <- read.table(file=lymphocyte.name,sep="\t",header=T)
slides.BRCA <- as.matrix(slides.BRCA)
################################################################################################
proportion_new <- slides.BRCA[8,]
cnames.BRCA <- colnames(BRCA)
cnames.proportion <- names(proportion_new)
if(nchar(cancername)==4){
  tumor.lable <- substr(cnames.BRCA,6,21)
}
if(nchar(cancername)==3){
  tumor.lable <- substr(cnames.BRCA,5,20)
}
if(nchar(cancername)==2){
  tumor.lable <- substr(cnames.BRCA,4,19)
}
proportion.lable <- substr(cnames.proportion,1,16)
match.location<-match(proportion.lable,tumor.lable)
na_match_idx<-which(is.na(match.location)==TRUE)
if(length(na_match_idx)>0){
  proportion_lable_tumor<-proportion.lable[-na_match_idx]
  proportion_afterfilter_matrix<-proportion_new[-na_match_idx]
}
if(length(na_match_idx)==0){
  proportion_lable_tumor<-proportion.lable
  proportion_afterfilter_matrix<-proportion_new
}
table<-table(proportion_lable_tumor)
table<-as.matrix(table)
table_duplicate_location<-which(table[,1]>1)
table_unduplicate_location<-which(table[,1]==1)
proportion_unduplicate_name<-names(table_unduplicate_location)
proportion_duplicate_name<-names(table_duplicate_location)
length(proportion_duplicate_name)
##################################################################################################
if(length(proportion_duplicate_name)==0){
  proportion_unduplicate_location<-match(proportion_unduplicate_name,proportion_lable_tumor)
  proportion_unduplicate_value<-matrix(proportion_afterfilter_matrix[proportion_unduplicate_location],
                                       ncol=1, nrow=length(proportion_unduplicate_name),
                                       dimnames=list(proportion_unduplicate_name,"value"))
  final_proportion_matrix<-proportion_unduplicate_value
}
###################################################################################################
if(length(proportion_duplicate_name)>1){
  final_sd<-c()
  for(i in (1:length(proportion_duplicate_name)))
  {
    temp_location<-which(proportion_lable_tumor==proportion_duplicate_name[i])
    temp_duplicate<-proportion_afterfilter_matrix[temp_location]
    temp_sd<-sd(temp_duplicate)
    final_sd<-c(final_sd,temp_sd)
  }
  hist(final_sd,freq=FALSE,breaks=100,col="red")
  axis(side = 1, at = c(1:40))
######################################################################################################  
  
  sd_delete<-which(final_sd>0.5)
  proportion_duplicate_name<-proportion_duplicate_name[-sd_delete]
  final_duplicate_value<-c()
  for(i in (1:length(proportion_duplicate_name)))
  {
    tmp_location<-which(proportion_lable_tumor==proportion_duplicate_name[i])
    tmp_duplicate<-proportion_afterfilter_matrix[tmp_location]
    tmp_mean<-mean(tmp_duplicate)
    final_duplicate_value<-c(final_duplicate_value,tmp_mean)
  }
  proportion_duplicate_value<-matrix(final_duplicate_value,ncol=1,nrow= length(final_duplicate_value),
                                     dimnames=list(proportion_duplicate_name,"value"))
  proportion_unduplicate_location<-match(proportion_unduplicate_name,proportion_lable_tumor)
  proportion_unduplicate_value<-matrix(proportion_afterfilter_matrix[proportion_unduplicate_location],
                                       ncol=1, nrow=length(proportion_unduplicate_name),
                                       dimnames=list(proportion_unduplicate_name,"value"))
  final_proportion_matrix<-rbind(proportion_duplicate_value,proportion_unduplicate_value)
}
#######################################################################################################
if(length(proportion_duplicate_name)==1){
  final_duplicate_value<-c()
  tmp_location<-which(proportion_lable_tumor==proportion_duplicate_name[1])
  tmp_duplicate<-proportion_afterfilter_matrix[tmp_location]
  tmp_mean<-mean(tmp_duplicate)
  final_duplicate_value<-c(final_duplicate_value,tmp_mean)
  proportion_duplicate_value<-matrix(final_duplicate_value,ncol=1,nrow= length(final_duplicate_value),
                                     dimnames=list(proportion_duplicate_name,"value"))
  proportion_unduplicate_location<-match(proportion_unduplicate_name,proportion_lable_tumor)
  proportion_unduplicate_value<-matrix(proportion_afterfilter_matrix[proportion_unduplicate_location],
                                       ncol=1, nrow=length(proportion_unduplicate_name),
                                       dimnames=list(proportion_unduplicate_name,"value"))
  final_proportion_matrix<-rbind(proportion_duplicate_value,proportion_unduplicate_value)
}
########################################################################################################
final_sample_name<-c(proportion_duplicate_name,proportion_unduplicate_name)
write.table(final_proportion_matrix,file=paste(cancername,"-yes",sep=""),sep="\t")
