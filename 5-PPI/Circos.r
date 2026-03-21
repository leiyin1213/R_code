########################################
# 1. 安装并加载R包
########################################
packages <- c("circlize", "biomaRt")

for(p in packages){
  if(!require(p, character.only = TRUE)){
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

########################################
# 2. 输入你的核心基因（hub genes）
########################################
genes <- read.table("Top4_intersection_MCCsorted.txt",sep="\t",header=T,check.names=F) [,1]

########################################
# 3. 获取基因染色体位置（Ensembl）
########################################
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_pos <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", 
                 "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

########################################
# 4. 数据整理
########################################
# 去掉非标准染色体（只保留1-22,X,Y）
gene_pos <- gene_pos[gene_pos$chromosome_name %in% c(1:22,"X","Y"), ]

# 转换为circlize格式
gene_df <- data.frame(
  chr = paste0("chr", gene_pos$chromosome_name),
  start = gene_pos$start_position,
  end = gene_pos$end_position,
  gene = gene_pos$hgnc_symbol
)

########################################
# 5. 绘制Circos图
########################################
pdf("circos_plot2.pdf", width=6, height=6)
circos.clear()
# ⭐ 关键：先设置参数
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "hg19")

# ===== 外圈颜色 =====
chr_colors <- rep(c("#4DBBD5","#E64B35","#00A087","#3C5488"), length.out=24)

circos.trackPlotRegion(
  ylim = c(0,1),
  track.height = 0.03,
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    chr_num = gsub("chr","", chr)
    
    idx = match(chr_num, c(1:22,"X","Y"))
    
    circos.rect(
      CELL_META$xlim[1], 0,
      CELL_META$xlim[2], 1,
      col = chr_colors[idx],
      border = NA
    )
  }
)

# ===== 基因标签（内侧）=====
circos.genomicLabels(
  gene_df,
  labels.column = 4,
  side = "inside",
  col = "black",
  line_col = "black",
  cex = 0.7,
  connection_height = mm_h(2)
)

dev.off()
 
