

expFile="singleGeneSurData.txt"      #表达输入文件
cliFile="clinical.txt"           #临床输入文件
methyFile="methyMatrix.txt"      #甲基化输入文件
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA数据库//2.甲基化数据//6.cliStat")      #设置工作目录

#读取表达文件，并对输入文件整理
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

#读取甲基化输入文件，并对数据进行整理
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#合并甲基化和表达数据
sameSample=intersect(row.names(exp), row.names(methy))
exp=exp[sameSample,]
methy=methy[sameSample,]
data=cbind(exp, methy=as.numeric(rowMeans(methy)))
data=data[,-2]
colnames(data)=c("expression","methylation")

#读取临床输入文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(cli),row.names(data))
cli=cli[sameSample,]
data=data[sameSample,]
rt=cbind(cli, data)

#数据分组
rt$expression=ifelse(rt$expression>median(rt$expression), "High", "Low")
rt$methylation=ifelse(rt$methylation>median(rt$methylation), "High", "Low")

#输出表达统计结果
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

#输出甲基化统计结果
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


