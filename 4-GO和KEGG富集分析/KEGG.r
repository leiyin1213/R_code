
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
BiocManager::install("ComplexHeatmap", update=FALSE, ask=FALSE)
BiocManager::install("cols4all")

devtools::install_github("davidsjoberg/ggsankey")
BiocManager::install("cowplot")
install.packages(c("shinyjs", "kableExtra", "colorblindcheck"))


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
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(cols4all)

library(tidyverse)
library(devtools)
library(ggsankey)
library(cowplot)
library(shinyjs)
library(kableExtra)
library(colorblindcheck)

setwd("D:/r/R_code/4-GO和KEGG富集分析")
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

kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

#保存结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GObarplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")


#GO圈图整理
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
#设置数目
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#画图
pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

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

# #绘制弦图
# # 修改后的代码段
# rt <- read.table("diff.Wilcoxon.txt", sep="\t", header=T, check.names=F) 
# KEGG <- read.table("KEGG.txt", sep="\t", header=T, check.names=F) 
# rownames(rt) <- rt$gene

# kegg1 <- strsplit(KEGG$geneID,"/")
# a1 <- data.frame(Term=KEGG$Description[1], Genes=kegg1[1])
# a2 <- data.frame(Term=KEGG$Description[2], Genes=kegg1[2])
# a3 <- data.frame(Term=KEGG$Description[3], Genes=kegg1[3])
# a4 <- data.frame(Term=KEGG$Description[4], Genes=kegg1[4])
# a5 <- data.frame(Term=KEGG$Description[5], Genes=kegg1[5])

# colnames(a1)[2] <- "gene"
# colnames(a2)[2] <- "gene"
# colnames(a3)[2] <- "gene"
# colnames(a4)[2] <- "gene"
# colnames(a5)[2] <- "gene"

# all1 <- rbind(a1,a2,a3,a4,a5)
# all1$value <- 1

# all2 <- spread(all1, key="Term", value="value", fill=0)
# rownames(all2) <- all2$gene
# all2 <- all2[-1]  # 移除第一列gene列

# # 确保基因顺序一致且维度匹配
# id <- intersect(rownames(all2), rownames(rt))
# print(paste("交集基因数量:", length(id)))

# # 只保留共同的基因
# all2 <- all2[id, ]
# rt_filtered <- rt[id, ]

# # 按logFC排序
# rt_filtered <- rt_filtered[order(rt_filtered$logFC, decreasing = T),]

# # 确保all2行顺序与排序后的rt一致
# all2 <- all2[rownames(rt_filtered), ]
# all2$logFC <- rt_filtered$logFC

# # 再次检查维度
# print(paste("最终all2行数:", nrow(all2)))
# print(paste("logFC向量长度:", length(all2$logFC)))
# pdf(file="KEGGGOChord.pdf", width=9, height=7)
# GOChord(all2,
#         space = 0.001,
#         gene.order = "logFC",
#         gene.size = 2,
#         gene.space = 0.25,
#         border.size = 0.1,
#         process.label = 7)
# dev.off()
# #Reactome 通路富集（ReactomePA）

# react_enrich <- enrichPathway(gene          = gene,
#                               organism      = "human",
#                               pvalueCutoff  = 1,
#                               qvalueCutoff  = 1,
#                               pAdjustMethod = "BH",
#                               readable      = TRUE)

# Reactome_res <- as.data.frame(react_enrich)
# Reactome_res_f <- Reactome_res %>% filter(pvalue < pvalueFilter & qvalue < qvalueFilter)

# write.table(Reactome_res_f, file="Reactome.txt", sep="\t", quote=FALSE, row.names=FALSE)

# showNum <- 20

# pdf("Reactome_barplot.pdf", width=10, height=7)
# print(barplot(react_enrich, showCategory=showNum, label_format=120, color=colorSel))
# dev.off()

# pdf("Reactome_dotplot.pdf", width=10, height=7)
# print(dotplot(react_enrich, showCategory=showNum, label_format=120, color=colorSel))
# dev.off()

# # Reactome通路之间的相似性网络（先计算pairwise相似度）
# react_sim <- pairwise_termsim(react_enrich)

# pdf("Reactome_emapplot.pdf", width=12, height=9)
# print(emapplot(react_sim, showCategory=30))
# dev.off()

# ==============================
#  绘制桑基图加气泡图
# ==============================
# 
#获取信息
kegg= KEGG[,c("Description","Count","pvalue","GeneRatio")]
#转化为分数转为小数
kegg$GeneRatio <- sapply(kegg$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#获取信息
sankey= KEGG[,c("Description","geneID")]
sankey <- sankey %>%
  separate_rows(geneID, sep = "/")

#数据处理
kegg2 <- kegg[length(rownames(kegg)):1,] 
kegg2 <- kegg2 %>%
  mutate(ymax = cumsum(Count)) %>%
  mutate(ymin = ymax -Count) %>%
  mutate(label = (ymin + ymax)/2)

#自定义主题
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))
#画图
p1 <- ggplot() +
  geom_point(data = kegg2,
             aes(x = -log10(pvalue),
                 y = label,
                 size = Count,
                 color = GeneRatio)) +
  scale_size_continuous(range=c(2,10)) +
  scale_y_continuous(expand = c(0,0.1),limits = c(0,52)) +
  scale_x_continuous(limits = c(0.1,ceiling(max(-log10(kegg2$pvalue)))+1)) +
  scale_colour_distiller(palette = "Reds", direction = 1) +
  labs(x = "-log10(Pvalue)",
       y = "") +
  theme_bw() +
  mytheme
# p1 <- dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel) +
#   theme(
#     axis.text.y = element_blank(),  # 去除y轴刻度标签
#     axis.ticks.y = element_blank()  # 去除y轴刻度线
#   )
p1

#数据处理
df <- sankey %>%
  make_long(geneID, Description)

#指定绘图顺序（转换为因子）：
df$node <- factor(df$node,levels = c(sankey$Description %>% unique()%>% rev(),
                                     sankey$geneID %>% unique() %>% rev()))

mycol <- c4a('rainbow_wh_rd',length(unique(df$node)))

#画图
p2 <- ggplot(df, aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80',
              node.fill = mycol, 
              smooth = 8,
              width = 0.08) +
  geom_sankey_text(size = 3.2,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')+ 
  theme(plot.margin = unit(c(0,8,0,0),units="cm"))
p2

#拼图(在p5的空白位置中插入p3)：
p3 <- ggdraw() + draw_plot(p2) + draw_plot(p1, scale = 0.5, x = 0.50, y=-0.17, width=0.68, height=1.3)
p3
ggsave("sankey_dot_plot.pdf", p3, width = 10, height = 7) 

