# ===== 0) 准备：feature genes（把你LASSO的结果填进来）=====
feature_genes <- c("LY96","TLR1","CYTIP","MS4A6A")

# ===== 1) 对齐样本顺序（很关键，不对齐会画错）=====
common_samples <- intersect(colnames(expr), names(group))
expr2 <- expr[, common_samples, drop = FALSE]
group2 <- group[common_samples]
group2 <- factor(group2, levels = c("Healthy","Disease"))  # 可按你需要调整顺序

# ===== 2) 提取4个基因的表达，转成长表 =====
library(dplyr)
library(tidyr)

df_long <- expr2[feature_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expr") %>%
  mutate(Group = group2[Sample])

# （可选但强烈推荐）如果不同平台/批次，先做标准化展示更美观：
# 每个基因做 z-score（只影响展示，不改变统计）
df_long <- df_long %>%
  group_by(Gene) %>%
  mutate(Expr_z = as.numeric(scale(Expr))) %>%
  ungroup()

# ===== 3) 统计检验：Wilcoxon（与原文一致）=====
stat_df <- df_long %>%
  group_by(Gene) %>%
  summarise(
    p = wilcox.test(Expr ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("p=", signif(p, 3)))

print(stat_df)

# ===== 4) 作图：箱线图 + 散点（对标 Fig.4D）=====
library(ggplot2)

ggplot(df_long, aes(x = Group, y = Expr_z, color = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = NULL, y = "Expression (z-score)")