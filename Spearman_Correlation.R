############################################################################################################

### Use the lymphocyte ratio to do spearman correlation.

############################################################################################################
tumormatrix_match<-match(final_sample_name,tumor.lable)
brca_matrix <- BRCA[,tumormatrix_match]
normalization<-function(x)
{
  
  temp_norm<-normalize.quantiles(x,copy=TRUE)
  rownames(temp_norm)<-rownames(x)
  colnames(temp_norm)<-colnames(x)
  return(temp_norm)
}
brca_norm <- normalization(brca_matrix)
rho<-c()
pvalue<-c()
each_gene_linear<-function(x)
{
  y<-final_proportion_matrix[,1]
  for(i in (1:nrow(x)))
  {
    vector_x<-x[i,]
    cor_result<-cor.test(vector_x,y,method = "spearman", conf.level = 0.95,alternative = "two.sided",exact = F)
    temp_rho<-cor_result$estimate
    temp_pvalue<-cor_result$p.value
    rho<-c(rho,temp_rho)
    pvalue<-c(pvalue,temp_pvalue)
  }
  resultset<-data.frame(rho,pvalue)
  colnames(resultset)<-c("rho","pvalue")
  rownames(resultset)<-rownames(x)
  return(resultset)
  
}
linear_brca_result <- each_gene_linear(brca_norm)
linear_brca_result <- as.matrix(linear_brca_result)
pvalue_signif_location<-which(linear_brca_result[,2]<0.05)
length(pvalue_signif_location)
pvalue_signif_brca<-linear_brca_result[pvalue_signif_location,]
top50 <- pvalue_signif_brca[order(pvalue_signif_brca[,2]),]
top50 <- top50[1:50,]
top10 <- top50[1:10,]
top500 <- pvalue_signif_brca[order(pvalue_signif_brca[,2]),]
top500 <- top500[1:500,]
write.table(top500,file=paste(cancername,"-top500",sep=""),sep="\t")
names.top50 <- rownames(top50)
names.top10 <- rownames(top10)
##############################################################################################################