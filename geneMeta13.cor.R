v

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(ggplot2)
library(ggpubr)
expFile="singleGeneSurData.txt"       #表达输入文件
methyFile="methyMatrix.txt"       #甲基化输入文件
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA数据库//2.甲基化数据//4.cor")       #设置工作目录

#读取表达输入文件，并对数据进行整理
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]
gene=colnames(exp)[1]

#读取甲基化输入文件，并对数据进行整理
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#定义相关性分析函数
bioCor=function(x=null, y=null, xlab=null, ylab=null, outFile=null){
	#相关性分析
	df1=as.data.frame(cbind(x,y))
	corT=cor.test(x,y,method="spearman")
	cor=corT$estimate
	pValue=corT$p.value
	p1=ggplot(df1, aes(x, y)) + 
				xlab(paste0(xlab," methylation"))+ylab(paste0(ylab," expression"))+
				geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				stat_cor(method = 'spearman', aes(x =x, y =y))
	#绘制相关性图形
	pdf(file=outFile, width=5, height=4.75)
	print(p1)
	dev.off()
}

#基因甲基化水平和基因表达的相关性分析
sameSample=intersect(row.names(exp), row.names(methy))
exp=exp[sameSample,]
methy=methy[sameSample,]
x=as.numeric(rowMeans(methy))
y=as.numeric(exp[,1])
bioCor(x=x, y=y, xlab=gene, ylab=gene, outFile="cor.gene.pdf")

#位点甲基化水平和基因表达的相关性分析
for(i in colnames(methy)){
	x=as.numeric(methy[,i])
	bioCor(x=x, y=y, xlab=i, ylab=gene, outFile=paste0("cor.", i, ".pdf"))
}

