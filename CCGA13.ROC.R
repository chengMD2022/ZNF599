

#install.packages("survivalROC")

library(survivalROC)
setwd("C://Users//Administrator//Desktop//������//ZNF599//1.TCGA���ݿ�//1.��������//13.ROC")      #���ù���Ŀ¼
rt=read.table("singleGeneSurData.txt",header=T,sep="\t",check.names=F,row.names=1)    #��ȡcox�ع�����ļ�
rocCol=c("red","green","blue")
aucText=c()

#����5���ROC����
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#����3���ROC����
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

#����1���ROC����
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()
