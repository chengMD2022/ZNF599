

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
methyFile="methyMatrix.txt"          #甲基化数据文件
cliFile="time.txt"     #生存输入文件
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA数据库//2.甲基化数据//5.methyOS")    #设置工作目录

#读取生存文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("OS.time","OS")]
colnames(cli)=c("futime","fustat")
cli=na.omit(cli)

#读取输入文件，提取位点的甲基化程度
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#合并数据
sameSample=intersect(row.names(cli), row.names(methy))
cli=cli[sameSample,]
methy=methy[sameSample,]
rt=cbind(cli, methy)
rt$futime=rt$futime/365

#生存分析
for(site in colnames(rt)[3:ncol(rt)]){
	group=ifelse(rt[,site]>median(rt[,site]), "High", "Low")
	diff=survdiff(Surv(futime, fustat) ~group, data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
			
	#绘制生存曲线
	fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
	surPlot=ggsurvplot(fit, 
			           data=rt,
			           conf.int=TRUE,
			           pval=pValue,
			           pval.size=6,
			           legend.labs=c("High", "Low"),
			           legend.title=site,
			           xlab="Time (years)",
			           ylab="Overall survival",
			           break.time.by = 1,
			           risk.table.title="",
			           palette=c("red", "blue"),
			           risk.table=T,
			           risk.table.height=.25)
	pdf(file=paste0("sur.",site,".pdf"), onefile=FALSE, width=8, height=5.5)
	print(surPlot)
	dev.off()
}

