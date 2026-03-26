# ============================================
# Figure 5C - ROC Curve
# ============================================

# ===== 1. 加载包 =====
library(tidyverse)
library(pROC)

# ===== 2. 读取数据（和你前面一样）=====
setwd("D:/r/R_code/6-lasso/GSE7158/单基因诊断")

expr <- read.table("GSE7158-normalize.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

disease <- read.table("Treat.txt")$V1
healthy <- read.table("Control.txt")$V1

# ===== 3. 分组 =====
group <- c(
  setNames(rep(1, length(disease)), disease),  # 1 = Disease
  setNames(rep(0, length(healthy)), healthy)   # 0 = Healthy
)

# ===== 4. 特征基因 =====
feature_genes <- read.table("intersectGenes.txt")$V1

# ===== 5. 对齐样本 =====
common_samples <- intersect(colnames(expr), names(group))

expr2 <- expr[feature_genes, common_samples]
group2 <- group[common_samples]

# ===== 6. 转置（样本为行）=====
data_mat <- t(expr2)
data_mat <- as.data.frame(data_mat)

# ===== 7. ROC分析（单基因）=====
library(ggplot2)

roc_list <- list()
auc_list <- c()

for (gene in feature_genes) {
  
  roc_obj <- roc(group2, data_mat[[gene]])
  
  roc_list[[gene]] <- roc_obj
  auc_list[gene] <- auc(roc_obj)
}

# ===== 8. 画图 =====
plot(NULL, xlim = c(1,0), ylim = c(0,1),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC Curve")

# ===== 颜色 =====
colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488")

i <- 1
for (gene in feature_genes) {
  
  lines(roc_list[[gene]],
        col = colors[i],
        lwd = 2)
  
  i <- i + 1
}

# ===== 图例（关键）=====
legend("bottomright",
       legend = paste0(names(auc_list),
                       " (AUC = ",
                       round(auc_list, 3), ")"),
       col = colors,
       lwd = 2,
       cex = 0.8)

# ===== 保存 =====
dev.copy(pdf, "Figure5C_ROC.pdf", width = 5, height = 5)
dev.off()