# ============================================
# Figure 5B - Violin Plot（论文复刻版）
# ============================================

# ===== 1. 加载包 =====
library(tidyverse)
library(ggplot2)
library(ggpubr)   # 🔥 关键包

# ===== 2. 设置路径 =====
setwd("D:/r/R_code/6-lasso/GSE7158/final")

expr_file <- "GSE7158-normalize.txt"
disease_file <- "Treat.txt"
healthy_file <- "Control.txt"

# ===== 3. 读取数据 =====
expr <- read.table(expr_file,
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

disease <- read.table(disease_file, stringsAsFactors = FALSE)$V1
healthy <- read.table(healthy_file, stringsAsFactors = FALSE)$V1

# ===== 4. 分组 =====
group <- c(
  setNames(rep("Disease", length(disease)), disease),
  setNames(rep("Healthy", length(healthy)), healthy)
)

# ===== 5. 特征基因 =====
feature_genes <- read.table("intersectGenes.txt",
                            stringsAsFactors = FALSE)$V1

# ===== 6. 对齐样本 =====
common_samples <- intersect(colnames(expr), names(group))

expr2 <- expr[, common_samples, drop = FALSE]
group2 <- factor(group[common_samples], levels = c("Healthy","Disease"))

# ===== 7. 转长表 =====
df_long <- expr2[feature_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expr") %>%
  mutate(Group = group2[Sample]) %>%
  filter(!is.na(Group))

# ===== 8. Z-score =====
df_long <- df_long %>%
  group_by(Gene) %>%
  mutate(Expr_z = as.numeric(scale(Expr))) %>%
  ungroup()

# ===== 9. 作图（🔥论文标准版）=====
for (g in unique(df_long$Gene)) {
  
  df_sub <- df_long %>% filter(Gene == g)
  
  # ===== 计算y轴上限（给括号留空间）=====
  y_max <- max(df_sub$Expr_z, na.rm = TRUE)
  
  p <- ggplot(df_sub, aes(x = Group, y = Expr_z, fill = Group)) +
    
    # 小提琴
    geom_violin(trim = FALSE, width = 0.7, alpha = 0.8, color = NA) +
    
    # 中位数线
    stat_summary(fun = median,
                 geom = "crossbar",
                 width = 0.3,
                 color = "black",
                 size = 0.6) +
    
    # 散点
    geom_jitter(width = 0.1,
                size = 1,
                color = "black",
                alpha = 0.8) +
    
    # 🔥 关键：括号 + p值（论文风格）
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",   # 显示 p = xxx
      label.y = y_max * 1.15,
      size = 4
    ) +
    
    scale_fill_manual(values = c("Healthy" = "#4DBBD5",
                                 "Disease" = "#E64B35")) +
    
    theme_classic() +
    theme(
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.background = element_rect(fill = "white", color = "black")
    ) +
    
    labs(
      title = g,
      x = NULL,
      y = "Expression (z-score)",
      fill = "Group"
    )
  
  print(p)
  
  ggsave(paste0("Figure5B_", g, ".pdf"),
         p,
         width = 10,
         height = 10)
}

