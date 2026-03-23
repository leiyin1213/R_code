# ===============================
# 1. 加载包
# ===============================
library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd("D:/r/R_code/6-lasso/GSE7158/final")
# ===============================
# 2. 读取分组信息（你上传的txt）
# ===============================
disease <- read.table("Treat.txt", stringsAsFactors = FALSE)$V1
healthy <- read.table("Control.txt", stringsAsFactors = FALSE)$V1

group <- data.frame(
  Sample = c(disease, healthy),
  Type = c(rep("Diseased", length(disease)),
           rep("Healthy", length(healthy)))
)

# ===============================
# 3. 读取表达矩阵
# ===============================
# ⚠️ 修改为你的文件名
expr <- read.table("GSE7158-normalize.txt", 
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t",
                   check.names = FALSE)

# ===============================
# 4. 提取4个基因
# ===============================
genes <-  read.table("intersectGenes.txt", stringsAsFactors = FALSE)$V1

# 检查基因是否存在
print(genes %in% rownames(expr))

# ===============================
# 5. 转换为长格式（作图关键）
# ===============================
data_long <- expr[genes, ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(group, by = "Sample")

# ===============================
# 6. 检查数据（避免踩坑）
# ===============================
# 是否有NA（说明样本没匹配上）
print(sum(is.na(data_long$Type)))

# ===============================
# 7. 绘制 Figure 4D
# ===============================
p <- ggplot(data_long, aes(x = Type, y = Expression, color = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Gene, scales = "free") +
  scale_color_manual(values = c("Diseased" = "#E64B35",
                                "Healthy" = "#4DBBD5")) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    text = element_text(size = 12)
  ) +
  labs(x = "", y = "Expression")

# 显示图
print(p)

# ===============================
# 8. 保存图片（论文用）
# ===============================
ggsave("Figure4D.pdf", p, width = 8, height = 6)