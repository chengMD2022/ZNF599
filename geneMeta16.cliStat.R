

expFile="singleGeneSurData.txt"      #���������ļ�
cliFile="clinical.txt"           #�ٴ������ļ�
methyFile="methyMatrix.txt"      #�׻��������ļ�
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA���ݿ�//2.�׻�������//6.cliStat")      #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

#��ȡ�׻��������ļ����������ݽ�������
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#�ϲ��׻����ͱ�������
sameSample=intersect(row.names(exp), row.names(methy))
exp=exp[sameSample,]
methy=methy[sameSample,]
data=cbind(exp, methy=as.numeric(rowMeans(methy)))
data=data[,-2]
colnames(data)=c("expression","methylation")

#��ȡ�ٴ������ļ�
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#�ϲ�����
sameSample=intersect(row.names(cli),row.names(data))
cli=cli[sameSample,]
data=data[sameSample,]
rt=cbind(cli, data)

#���ݷ���
rt$expression=ifelse(rt$expression>median(rt$expression), "High", "Low")
rt$methylation=ifelse(rt$methylation>median(rt$methylation), "High", "Low")

#�������ͳ�ƽ��
cliStatOut=data.frame()
for(i in 1:ncol(rt)){
	nameStat=colnames(rt)[i]
	tableStat=table(rt[,c(nameStat,"expression")])
	tableStatSum=cbind(Total=rowSums(tableStat),tableStat)
	tableStatRatio=prop.table(tableStatSum,2)
	tableStatRatio=round(tableStatRatio*100,2)
	tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
	tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
	pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
	pValueStat=round(pStat$p.value,4)
	pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
	tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
	cliStatOut=rbind(cliStatOut,tableStatOut)
}
write.table(cliStatOut,file="cliStat.exp.xls",sep="\t",quote=F,row.names=F)

#����׻���ͳ�ƽ��
cliStatOut=data.frame()
for(i in 1:ncol(rt)){
	nameStat=colnames(rt)[i]
	tableStat=table(rt[,c(nameStat,"methylation")])
	tableStatSum=cbind(Total=rowSums(tableStat),tableStat)
	tableStatRatio=prop.table(tableStatSum,2)
	tableStatRatio=round(tableStatRatio*100,2)
	tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
	tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
	pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
	pValueStat=round(pStat$p.value,4)
	pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
	tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
	cliStatOut=rbind(cliStatOut,tableStatOut)
}
write.table(cliStatOut,file="cliStat.methy.xls",sep="\t",quote=F,row.names=F)

