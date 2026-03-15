
setwd("D:/r/实验数据/5-PPI")
library(readr)
# 假设 MCC 和 Degree 的结果分别保存在两个数据框中

mcc <- read_csv("mcc715.csv", )
deg <- read_csv("Degree715.csv")

# 看看列名长什么样（先跑一次）
colnames(mcc)
colnames(deg)

pick_gene_col <- function(df){
  candidates <- c("Gene","gene","Name","name","Node","node","Protein","protein","ID","id","node_name")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  # 如果都没命中，就默认第一列是基因名
  return(colnames(df)[1])
}

gene_col_mcc <- pick_gene_col(mcc)
gene_col_deg <- pick_gene_col(deg)

mcc_genes <- unique(trimws(mcc[[1]]))
deg_genes <- unique(trimws(deg[[1]]))

common_genes <- intersect(mcc_genes, deg_genes)
common_genes
length(common_genes)

pick_score_col <- function(df){
  candidates <- c("Score","score","MCC","Degree","rank","Rank")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  return(NULL)
}

score_col_mcc <- pick_score_col(mcc)

mcc_common <- mcc[trimws(mcc[[gene_col_mcc]]) %in% common_genes, ]

if (!is.null(score_col_mcc)) {
  mcc_common[[score_col_mcc]] <- as.numeric(mcc_common[[score_col_mcc]])
  mcc_common <- mcc_common[order(mcc_common[[score_col_mcc]], decreasing = TRUE), ]
}

top4 <- head(trimws(mcc_common[[gene_col_mcc]]), 10)
write.csv(data.frame(Gene=mcc_genes), "Top4_intersection_MCCsorted.csv",
          row.names = FALSE, quote = FALSE)
write.table(top4, "Top4_intersection_MCCsorted.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

