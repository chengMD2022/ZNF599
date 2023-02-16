

#install.packages("plyr")
#install.packages("ggpubr")


#引用包
library(plyr)
library(reshape2)
library(ggpubr)
methyFile="methyMatrix.txt"     #甲基化输入文件
outFile="boxplot.pdf"           #输出图片文件
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA数据库//2.甲基化数据//3.boxplot")       #设置工作目录
methy=read.table(methyFile, header=T, sep="\t", check.names=F,row.names=1)    #读取输入文件
methy=t(methy)

#把数据转换成ggplot2输入文件
rt=melt(methy)
colnames(rt)=c("Id","Site","Expression")

#定义输出图片的排序方式
med=ddply(rt,"Site",summarise,med=median(Expression))
rt$Site=factor(rt$Site, levels=med[order(med[,"med"],decreasing = T), "Site"])

#绘制箱线图
p=ggboxplot(rt, x="Site", y="Expression", color="Site",
       palette = rainbow(ncol(methy)),
       ylab="Methylation level (beta value)",
       xlab="Methylation site",
       width=0.6,
       #add = "jitter",                  #绘制每个样品的散点
       legend = "right")
pdf(file=outFile, width=8, height=6)    #输出图片文件
p+rotate_x_text(45)
dev.off()

