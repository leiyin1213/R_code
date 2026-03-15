
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")
install.packages("ggplot2")
install.packages("stringi")
install.packages("GOplot")
BiocManager::install("ReactomePA", update=FALSE, ask=FALSE)

library(R.utils)

#加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	
library(GOplot)
library(tidyr)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(ReactomePA)

setwd("D:/r/实验数据/4-GO和KEGG富集分析")
#读入
input_diff =read.table("intersectGenes1.txt",sep="\t",header=T,check.names=F) 
input_gene <- input_diff[,1]
#去除重复基因
input_gene=unique(as.vector(input_gene))

#将gene symbol转换为基因id
#org.Hs.eg为人的物种
#https://www.jianshu.com/p/84e70566a6c6
entrezIDs=mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
#去除基因id为NA的基因
gene=entrezIDs[entrezIDs!="NA"]
#去除多个ID
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#筛选条件
pvalueFilter=0.05
qvalueFilter=1       

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

#GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

#保存
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#显示的数目
showNum=20

#柱状图
pdf(file="GObarplot.pdf", width=20, height=15)
# bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
bar_plot <- barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel)
bar <- bar_plot + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.pdf", width=20, height=15)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

barplot(kk)  #富集柱形图
dotplot(kk)  #富集气泡图
cnetplot(kk) #网络图展示富集功能和基因的包含关系
emapplot(kk) #网络图展示各富集功能之间共有基因关系
heatplot(kk) #热图展示富集功能和基因的包含关系

#筛选条件
pvalueFilter=0.05
qvalueFilter=1       

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

#KEGG
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(input_gene[match(strsplit(x,"/")[[1]],as.character(entrezIDs))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

#保存
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=20

#柱状图
pdf(file="KEGGbarplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

#绘制弦图
# 修改后的代码段
rt <- read.table("diff.Wilcoxon.txt", sep="\t", header=T, check.names=F) 
KEGG <- read.table("KEGG.txt", sep="\t", header=T, check.names=F) 
rownames(rt) <- rt$gene

kegg1 <- strsplit(KEGG$geneID,"/")
a1 <- data.frame(Term=KEGG$Description[1], Genes=kegg1[1])
a2 <- data.frame(Term=KEGG$Description[2], Genes=kegg1[2])
a3 <- data.frame(Term=KEGG$Description[3], Genes=kegg1[3])
a4 <- data.frame(Term=KEGG$Description[4], Genes=kegg1[4])
a5 <- data.frame(Term=KEGG$Description[5], Genes=kegg1[5])

colnames(a1)[2] <- "gene"
colnames(a2)[2] <- "gene"
colnames(a3)[2] <- "gene"
colnames(a4)[2] <- "gene"
colnames(a5)[2] <- "gene"

all1 <- rbind(a1,a2,a3,a4,a5)
all1$value <- 1

all2 <- spread(all1, key="Term", value="value", fill=0)
rownames(all2) <- all2$gene
all2 <- all2[-1]  # 移除第一列gene列

# 确保基因顺序一致且维度匹配
id <- intersect(rownames(all2), rownames(rt))
print(paste("交集基因数量:", length(id)))

# 只保留共同的基因
all2 <- all2[id, ]
rt_filtered <- rt[id, ]

# 按logFC排序
rt_filtered <- rt_filtered[order(rt_filtered$logFC, decreasing = T),]

# 确保all2行顺序与排序后的rt一致
all2 <- all2[rownames(rt_filtered), ]
all2$logFC <- rt_filtered$logFC

# 再次检查维度
print(paste("最终all2行数:", nrow(all2)))
print(paste("logFC向量长度:", length(all2$logFC)))
pdf(file="KEGGGOChord.pdf", width=9, height=7)
GOChord(all2,
        space = 0.001,
        gene.order = "logFC",
        gene.size = 2,
        gene.space = 0.25,
        border.size = 0.1,
        process.label = 7)
dev.off()
#Reactome 通路富集（ReactomePA）

react_enrich <- enrichPathway(gene          = gene,
                              organism      = "human",
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              readable      = TRUE)

Reactome_res <- as.data.frame(react_enrich)
Reactome_res_f <- Reactome_res %>% filter(pvalue < pvalueFilter & qvalue < qvalueFilter)

write.table(Reactome_res_f, file="Reactome.txt", sep="\t", quote=FALSE, row.names=FALSE)

showNum <- 20

pdf("Reactome_barplot.pdf", width=10, height=7)
print(barplot(react_enrich, showCategory=showNum, label_format=120, color=colorSel))
dev.off()

pdf("Reactome_dotplot.pdf", width=10, height=7)
print(dotplot(react_enrich, showCategory=showNum, label_format=120, color=colorSel))
dev.off()

# Reactome通路之间的相似性网络（先计算pairwise相似度）
react_sim <- pairwise_termsim(react_enrich)

pdf("Reactome_emapplot.pdf", width=12, height=9)
print(emapplot(react_sim, showCategory=30))
dev.off()
