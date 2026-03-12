
#install.packages("pheatmap")

#install.packages("stringr")

#install.packages("ggplot2")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("devtools")
#devtools::install_github("BioSenior/ggVolcano")

library(limma)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggVolcano)

#设置工作目录
setwd("D:/r/实验数据/2-Degs/GSE16134")

#读取输入文件
data=read.table("GSE16134.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
#data=data[rowMeans(data)>1,]
boxplot(data.frame(data),col="#4DBBD5")
#标准化
data=normalizeBetweenArrays(data)
boxplot(data.frame(data),col="#4DBBD5")
dev.off()

#导出
write.table(data.frame(ID=rownames(data),data),file="normalize.txt", sep="\t", quote=F, row.names = F)

#读入Control样本
Control=read.table("Control.txt", header=F, sep="\t", check.names=F)

#读入Treat样本
Treat=read.table("Treat.txt", header=F, sep="\t", check.names=F)

#设置分组信息
conNum=length(rownames(Control))
treatNum=length(rownames(Treat))
Type=c(rep(1,conNum), rep(2,treatNum))

#按照正常，治疗排序
data1 = data[,Control[,1]]
data2 = data[,Treat[,1]]
data = cbind(data1,data2)

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  #整理时，已经是log2后的
#  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  logFC=treatGeneMeans-conGeneMeans
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
#  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
 # }
}

pValue=outTab[,"pValue"]
#fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
fdr = p.adjust(pValue, method="fdr")
outTab=cbind(outTab, fdr=fdr)

write.table(outTab,file="all.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

logFCfilter<- log2(1)
fdrFilter=0.05
#输出差异基因
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & 
                  as.numeric(as.vector(outTab$pValue))<fdrFilter),]
# outDiff = outTab[ abs(logFC)>logFCfilter & pValue<0.05, ]
write.table(outDiff,file="diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

#热图
geneNum=50     
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=6.5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("#4DBBD5",5), "white", rep("#E64B35",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#火山图
pdf(file="vol.pdf", width=5, height=5)
xMax=2.5
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(-3,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1.2)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#E64B35",cex=1.5)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#4DBBD5",cex=1.5)
abline (h=-log10(0.05), lty=2, col="grey70")
abline (v=c(-logFCfilter, logFCfilter), lty=2, col="grey70")
abline(v=0, lty=2, col="grey40")
dev.off()

