

#install.packages("plyr")
#install.packages("ggpubr")


#���ð�
library(plyr)
library(reshape2)
library(ggpubr)
methyFile="methyMatrix.txt"     #�׻��������ļ�
outFile="boxplot.pdf"           #���ͼƬ�ļ�
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA���ݿ�//2.�׻�������//3.boxplot")       #���ù���Ŀ¼
methy=read.table(methyFile, header=T, sep="\t", check.names=F,row.names=1)    #��ȡ�����ļ�
methy=t(methy)

#������ת����ggplot2�����ļ�
rt=melt(methy)
colnames(rt)=c("Id","Site","Expression")

#�������ͼƬ������ʽ
med=ddply(rt,"Site",summarise,med=median(Expression))
rt$Site=factor(rt$Site, levels=med[order(med[,"med"],decreasing = T), "Site"])

#��������ͼ
p=ggboxplot(rt, x="Site", y="Expression", color="Site",
       palette = rainbow(ncol(methy)),
       ylab="Methylation level (beta value)",
       xlab="Methylation site",
       width=0.6,
       #add = "jitter",                  #����ÿ����Ʒ��ɢ��
       legend = "right")
pdf(file=outFile, width=8, height=6)    #���ͼƬ�ļ�
p+rotate_x_text(45)
dev.off()
