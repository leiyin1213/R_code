BiocManager::install("DOSE")
BiocManager::install("pathview")
install.packages("clusterProfiler")



library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("DOSE")

setwd("D:/r/实验数据/4-GO和KEGG富集分析")

pvalueFilter=0.05        
qvalueFilter=1       
showNum=25

rt=read.table("intersectGenes.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichDO(gene =gene, ont =  "HDO",pvalueCutoff = 1, qvalueCutoff = 1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="DO.xls",sep="\t",quote=F,row.names = F)
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}
pdf(file="DO_barplot.pdf",width =9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()
pdf(file="DO_bubble.pdf",width =9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

# 提取数据并转换 GeneRatio 为数值
df <- as.data.frame(kk)

# 将 GeneRatio 字符串 "a/b" 转换为数值 a/b
df$GeneRatio <- sapply(df$GeneRatio, function(x) {
  if (is.na(x)) return(NA_real_)
  parts <- strsplit(as.character(x), "/")[[1]]
  if (length(parts) == 2) {
    return(as.numeric(parts[1]) / as.numeric(parts[2]))
  } else {
    return(NA_real_)
  }
})

# 将 pvalue 转换为数值
df$pvalue <- as.numeric(as.character(df$pvalue))

# 计算 GenePercent
df$GenePercent <- df$GeneRatio * 100

# 检查是否成功
print(head(df[, c("Description", "GeneRatio", "GenePercent", "pvalue")]))
# 输出 PDF
pdf(file = "DO_barplot2.pdf", width = 9, height = 7)
# 创建水平条形图
p <- ggplot(df, aes(x = GenePercent, y = Description, fill = pvalue)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_gradient(low = "steelblue", high = "darkblue",
                      limits = c(min(df$pvalue, na.rm = TRUE), max(df$pvalue, na.rm = TRUE))) +
  labs(x = "Gene Percent(%)", y = "Doterm", fill = "pvalue") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 50)) +
  geom_text(aes(label = paste(round(GenePercent, 1), " (", format(pvalue, scientific = TRUE), ")", sep = "")),
            size = 3, hjust = -0.2, vjust = 0.5) +
  coord_flip()

# 输出 PDF
#pdf(file = "DO_barplot2.pdf", width = 9, height = 7)
print(p)
dev.off()

