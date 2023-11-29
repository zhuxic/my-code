

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #���ð�
setwd("D:\\melanoma\\anoikis\\48.preMaftools")      #���ù���Ŀ¼

#��ȡ�����ļ�����Ϊ�ٲ�ͼ��ע����Ϣ
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#��ȡ����ͻ����ļ�
geneNum=20     #�����ٲ�ͼ����չʾ�������Ŀ
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#����ע�͵���ɫ
ann_colors=list()
col=c("blue", "red")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#���Ƶͷ������ٲ�ͼ
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#���Ƹ߷������ٲ�ͼ
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


