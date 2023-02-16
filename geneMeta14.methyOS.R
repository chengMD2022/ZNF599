

#install.packages("survival")
#install.packages("survminer")


#���ð�
library(survival)
library(survminer)
methyFile="methyMatrix.txt"          #�׻��������ļ�
cliFile="time.txt"     #���������ļ�
setwd("C://Users//Administrator//Desktop//ZNF599//1.TCGA���ݿ�//2.�׻�������//5.methyOS")    #���ù���Ŀ¼

#��ȡ�����ļ�
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("OS.time","OS")]
colnames(cli)=c("futime","fustat")
cli=na.omit(cli)

#��ȡ�����ļ�����ȡλ��ļ׻����̶�
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy)

#�ϲ�����
sameSample=intersect(row.names(cli), row.names(methy))
cli=cli[sameSample,]
methy=methy[sameSample,]
rt=cbind(cli, methy)
rt$futime=rt$futime/365

#�������
for(site in colnames(rt)[3:ncol(rt)]){
	group=ifelse(rt[,site]>median(rt[,site]), "High", "Low")
	diff=survdiff(Surv(futime, fustat) ~group, data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
			
	#������������
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
