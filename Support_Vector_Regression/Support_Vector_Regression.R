#############################################################################################################

### Do the SVR

#############################################################################################################
final_proportion_matrix <- log2(final_proportion_matrix)
index <- match(names.top50,rownames(BRCA))
subdata <- BRCA[index,]
samples.name <- colnames(subdata)
subdata <- normalize.quantiles(subdata,copy = TRUE)
subdata <- t(subdata)
subdata <- scale(subdata)
#############################################################################################################
sample.index <- match(final_sample_name,tumor.lable)
samples.slides.yes <- samples.name[sample.index]
subdata.1 <- t(subdata)[,sample.index]
subdata.1 <- t(subdata.1)
na.index <- which(is.na(match(substr(samples.name,6,21),final_sample_name))==TRUE)
samples.slides.no <- samples.name[na.index]
subdata.2 <- subdata[na.index,]
rmse <- function(error)
{
  sqrt(mean(error^2))
}
library(e1071)
lable <- character()
for(i in 1:50){
  lable[i] <- paste("x",i,sep="")
}
subdata.1 <- cbind(final_proportion_matrix,subdata.1)
subdata.1 <- as.data.frame(subdata.1)
colnames(subdata.1)[2:51] <- lable
model <- svm(value ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+
               x23+x24+x25+x26+x27+x28+x29+x30+x31+x32+x33+x34+x35+x36+x37+x38+x39+x40+x41+x42+x43+
               x44+x45+x46+x47+x48+x49+x50, data = subdata.1)
colnames(subdata.2) <- lable
subdata.2 <- as.data.frame(subdata.2)
predictedY <- predict(model, newdata = subdata.2)
##############################################################################################################
tuneResult <- tune(svm, value ~ ., data = subdata.1, ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
tunedModel <- tuneResult$best.model
tunedModelY <- predict(tunedModel, newdata = subdata.2)
write.table(final_proportion_matrix,paste(cancername,"-final.proportion.matrix",sep=""),sep="\t")
write.table(samples.slides.yes,paste(cancername,"-samples.slides.yes",sep=""),sep="\t")
write.table(tunedModelY,paste(cancername,"-tunedModelY",sep=""),sep="\t")
write.table(samples.slides.no,paste(cancername,"-samples.slides.no",sep=""),sep="\t")
write.table(names.top50,file=paste(cancername,"-top50",sep=""),sep="\t")

a <- as.vector(predict(tunedModel, newdata = subdata.1))
b <- as.vector(final_proportion_matrix)
compare <- data.frame(True = b, Predict = a)
write.table(compare,file=paste(cancername,"_compare",sep=""),sep="\t")
png(file=paste(cancername,".png",sep=""))
plot(a,b,pch=19,xlab = "Predicted",ylab = "True")
dev.off()
cor.test(a,b,method="pearson",conf.level=0.95,alternative="two.sided")
dim(subdata)
dim(subdata.1)
dim(subdata.2)
