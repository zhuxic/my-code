

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")           #���ð�
inputFile="merge.txt"      #���������ļ�
setwd("D:\\melanoma\\anoikis\\44.CIBERSORT")      #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(inputFile, header=T, sep="\t", check.names=F)        #��ȡ�����ļ�
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#��������õ�����
out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #����ļ�

#����CIBERSORT���õ�����ϸ���������
source("prgTME44.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)


