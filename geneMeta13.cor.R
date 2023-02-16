v

#install.packages("ggplot2")
#install.packages("ggpubr")


#���ð�
library(ggplot2)
library(ggpubr)
expFile="singleGeneSurData.txt"       #���������ļ�
methyFile="methyMatrix.txt"       #�׻��������ļ�
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA���ݿ�//2.�׻�������//4.cor")       #���ù���Ŀ¼

#��ȡ���������ļ����������ݽ�������
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]
gene=colnames(exp)[1]

#��ȡ�׻��������ļ����������ݽ�������
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#��������Է�������
bioCor=function(x=null, y=null, xlab=null, ylab=null, outFile=null){
	#����Է���
	df1=as.data.frame(cbind(x,y))
	corT=cor.test(x,y,method="spearman")
	cor=corT$estimate
	pValue=corT$p.value
	p1=ggplot(df1, aes(x, y)) + 
				xlab(paste0(xlab," methylation"))+ylab(paste0(ylab," expression"))+
				geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				stat_cor(method = 'spearman', aes(x =x, y =y))
	#���������ͼ��
	pdf(file=outFile, width=5, height=4.75)
	print(p1)
	dev.off()
}

#����׻���ˮƽ�ͻ�����������Է���
sameSample=intersect(row.names(exp), row.names(methy))
exp=exp[sameSample,]
methy=methy[sameSample,]
x=as.numeric(rowMeans(methy))
y=as.numeric(exp[,1])
bioCor(x=x, y=y, xlab=gene, ylab=gene, outFile="cor.gene.pdf")

#λ��׻���ˮƽ�ͻ�����������Է���
for(i in colnames(methy)){
	x=as.numeric(methy[,i])
	bioCor(x=x, y=y, xlab=i, ylab=gene, outFile=paste0("cor.", i, ".pdf"))
}
