# ============================================
# Figure 5C - Correlation Plot
# ============================================

# ===== 1. 加载包 =====
library(tidyverse)
library(ggpubr)

# ===== 2. 读取表达矩阵 =====
setwd("D:/r/R_code/6-lasso/GSE7158/final")

expr <- read.table("GSE7158-normalize.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

# ===== 3. 选择两个基因（按论文改）=====
gene_x <- "IFI35"
gene_y <- "CXCL8"

# ===== 4. 提取表达 =====
df <- data.frame(
  x = as.numeric(expr[gene_x, ]),
  y = as.numeric(expr[gene_y, ])
)

# ===== 5. 计算相关性 =====
cor_test <- cor.test(df$x, df$y, method = "pearson")

r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)

label_text <- paste0("R = ", r_val, ", p = ", p_val)

# ===== 6. 作图 =====
p <- ggplot(df, aes(x = x, y = y)) +
  
  # 🔥 散点
  geom_point(color = "black", size = 2, alpha = 0.8) +
  
  # 🔥 回归线 + 置信区间
  geom_smooth(method = "lm",
              color = "red",
              fill = "red",
              alpha = 0.2,
              size = 1.2) +
  
  # 🔥 左上角标注（关键！）
  annotate("text",
           x = min(df$x),
           y = max(df$y),
           label = label_text,
           hjust = 0,
           size = 5) +
  
  # 🔥 论文风格
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  
  labs(
    x = gene_x,
    y = gene_y,
    title = "PD-GSE10334"
  )

# ===== 7. 输出 =====
print(p)

ggsave("Figure5C_correlation.pdf", p, width = 5, height = 5)

