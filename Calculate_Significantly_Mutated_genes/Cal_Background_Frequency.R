a <- read.csv(file = "motiffreqtypes.csv", header = T)
col1<-as.vector(a[,1])
motif1<-strsplit(col1," ")
snp<-c()
conten<-c()
for(i in 1:96){
  a1<-motif1[[i]][1]
  a2<-motif1[[i]][2]
  snp<-c(snp,a1)
  conten<-c(conten,a2)
}
snp[which(snp=="CA")]="CAGT"
snp[which(snp=="CG")]="CGGC"
snp[which(snp=="CT")]="CTGA"
snp[which(snp=="TA")]="TAAT"
snp[which(snp=="TC")]="TCAG"
snp[which(snp=="TG")]="TGAC"
motif<-data.frame(snp,conten)
infile<-read.csv("TCGA.SKCM.varscan.txt",header = T, sep = '\t')
s_n<-length(unique(infile$Tumor_Sample_Barcode))
indata<-infile[,c(1,5,6,7,9,10,11,12,13,16,55,56,57,112)]
indata<-indata[which(indata$Variant_Type=="SNP"),]
tri_motif<-substr(indata$CONTEXT,5,7)
substr(tri_motif,2,2)<-"."
subsit<-paste0(indata$Reference_Allele,indata$Tumor_Seq_Allele2)
num<-c()
for (i in 1:96){
  ind1<-which(subsit==substr(motif[i,1],1,2))
  ind2<-which(subsit==substr(motif[i,1],3,4))
  subsit_match<-c(ind1,ind2)
  conten_match<-which(tri_motif==as.character(motif[i,2]))
  match_index<-intersect(subsit_match,conten_match)
  num<-c(num,length(match_index))
}
num<-num/sum(num)
out<-data.frame(snv,conten,num)
skcm <- c(out$num,out$num)
back <- read.table(file = "background_cancer.txt", sep = '\t', header = T)
back$SKCM <- skcm
write.table(back,file="final_background.txt",sep="\t")
