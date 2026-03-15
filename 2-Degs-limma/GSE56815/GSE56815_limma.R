# 引用包
library(limma)
library(pheatmap)

gse <- "GSE56815"
base_path <- "D:/r/实验数据/2-Degs-limma/"
full_path <- file.path(base_path, gse)
setwd(full_path)

# 读取输入文件（已处理 series_matrix 导出的表达矩阵）
type <- ".txt"
expr <- read.table(paste0(gse, type), header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)

# 转为数值矩阵
expr <- as.matrix(expr)
mode(expr) <- "numeric"

# ========== 关键：检查是否需要 log2 ==========
qx <- quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
need_log2 <- (qx[6] > 100) || (qx[5] > 50 && qx[2] > 0)  # 经验规则
if (need_log2) {
  expr <- log2(expr + 1)
  message("检测到表达值偏大，已执行 log2(x+1)。")
} else {
  message("检测到表达值已在 log2 范围，跳过 log2。")
}

# ========== 可选：同平台芯片进一步分位数归一化 ==========
# series_matrix 可能已归一化；若你确定已归一化，可注释掉下一行
expr <- normalizeBetweenArrays(expr, method = "quantile")

# 导出归一化矩阵
write.table(data.frame(ID = rownames(expr), expr),
            file = paste0(gse, "-normalize.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 读入样本列表（每行一个样本名）
Control <- read.table("Control.txt", header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
Treat <- read.table("Treat.txt", header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)

control_samples <- Control[, 1]
treat_samples <- Treat[, 1]

# ========== 样本名检查 ==========
stopifnot(all(control_samples %in% colnames(expr)))
stopifnot(all(treat_samples %in% colnames(expr)))

# 提取并合并（按 Control 后 Treat 的顺序）
expr_sub <- expr[, c(control_samples, treat_samples), drop = FALSE]
conNum <- length(control_samples)
treatNum <- length(treat_samples)

# ========== limma 差异分析 ==========
group <- factor(c(rep("Control", conNum), rep("Treat", treatNum)),
                levels = c("Control", "Treat"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_sub, design)
cont.matrix <- makeContrasts(Treat - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ========== 导出全基因结果：用 BH（FDR） ==========
allDiff <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")

# ========== 新增：提取并添加SE值到结果表 ==========
# 取出当前对比（系数）对应的未缩放标准误
se_unscaled <- fit2$stdev.unscaled[, 1]   # 1 是第1个contrast/coef，按你实际对比列调整

# 用后验方差（与 moderated t 一致）
se_values <- se_unscaled * sqrt(fit2$s2.post)

# 给 se_values 加上基因名，保证能match
names(se_values) <- rownames(fit2$stdev.unscaled)

# 合并到 allDiff（allDiff 的行名必须是基因名）
allDiff$SE <- se_values[rownames(allDiff)]

write.table(allDiff, file = paste0(gse, "-all.gene.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# ========== 筛选阈值 ==========
logFCfilter <- log2(1)  # 你原来是 log2(1)=0，这等于不筛FC；这里建议至少 1 或 0.5
adjPfilter <- 0.05

# 如果你的数据确实很弱，可以先用 adjPfilter=0.1 或改用 P.Value 做探索
diffSig <- subset(allDiff, abs(logFC) > logFCfilter & P.Value < adjPfilter)
write.table(diffSig, file = paste0(gse, "-diff.gene.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# 区分上调和下调基因
upRegulated <- subset(allDiff, logFC > logFCfilter & P.Value < adjPfilter)
downRegulated <- subset(allDiff, logFC < -logFCfilter & P.Value < adjPfilter)

# 输出差异基因表达矩阵
diffGeneExp <- expr_sub[rownames(diffSig), , drop = FALSE]
write.table(data.frame(ID = rownames(diffGeneExp), diffGeneExp),
            file = paste0(gse, "-diffGeneExp.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 分别输出上调和下调基因名
write.table(rownames(upRegulated),
            file = paste0(gse, "-upGenename.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(downRegulated),
            file = paste0(gse, "-downGenename.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 也可以输出总的差异基因名（保持原有功能）
write.table(rownames(diffSig),
            file = paste0(gse, "-diffGenename.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========== 热图 ==========

diffSig2 <- diffSig[order(diffSig$logFC), , drop = FALSE]
diffGeneName <- rownames(diffSig2)
diffLength <- length(diffGeneName)
geneNum <-length(diffGeneName)
if (diffLength == 0) {
  message("没有筛到差异基因（FDR阈值过严或信号弱）。你可以先把 adjPfilter 放宽到 0.1 或改用 P.Value 做探索。")
} else {
  hmGene <- if (diffLength > (2 * geneNum)) {
    diffGeneName[c(1:geneNum, (diffLength - geneNum + 1):diffLength)]
  } else diffGeneName
  
  hmExp <- expr_sub[hmGene, , drop = FALSE]
  
  ann <- data.frame(Group = group)
  rownames(ann) <- colnames(expr_sub)
  
  # 分组颜色
  ann_colors <- list(
    Group = c(Control = "#18C3D6", Treat = "#F29CA3")
  )
  
  # 更接近参考图的蓝白红配色
  my_color <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # 固定色阶范围
  my_breaks <- seq(-2, 2, length.out = 101)
  
  pdf(file = paste0(gse, "-heatmap_compact.pdf"), width = 7.2, height = 4.2)
  
  # 你的目标基因
  target.genes <- c("DAPK2", "MAP1LC3B", "TUBA8", "RAB2A")
  
  # 假设表达矩阵叫 heatmap_matrix
  # 行 = 基因，列 = 样本
  # 先检查目标基因是否存在
  target.genes <- intersect(target.genes, rownames(hmExp))
  
  # 生成行标签：只有目标基因显示名字，其他为空
  lab_row <- ifelse(rownames(hmExp) %in% target.genes,
                    rownames(hmExp), "")
  pheatmap(hmExp,
           scale = "row",
           color = my_color,
           breaks = my_breaks,
           annotation_col = ann,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,   # 关键：隐藏基因名
           show_colnames = FALSE,   # 不显示样本名
           border_color = NA, 
           labels_row = lab_row,# 去边框更干净
           fontsize = 8,
           fontsize_row = 10,
           fontsize_col = 35)  
  
  dev.off()
}

# ========== 火山图 ==========
# 用 P.Value 画更直观；想用 FDR 就把 P.Value 换成 adj.P.Val
pdf(file = paste0(gse, "-vol.pdf"), width = 5, height = 5)

xMax <- max(2, ceiling(max(abs(allDiff$logFC), na.rm = TRUE)))
yVal <- -log10(allDiff$P.Value)
yMax <- max(yVal[is.finite(yVal)], na.rm = TRUE) + 1

plot(allDiff$logFC, yVal,
     xlab = "log2(Fold Change)", ylab = "-log10(P.Value)",
     main = paste0("Volcano: ", gse),
     ylim = c(0, yMax), xlim = c(-xMax, xMax),
     pch = 20, cex = 1.0)

# 标记显著点（这里用 FDR）
sigUp <- subset(allDiff, adj.P.Val < adjPfilter & logFC > logFCfilter)
sigDn <- subset(allDiff, adj.P.Val < adjPfilter & logFC < -logFCfilter)
points(sigUp$logFC, -log10(sigUp$P.Value), pch = 20, col = "#E64B35", cex = 1.2)
points(sigDn$logFC, -log10(sigDn$P.Value), pch = 20, col = "#4DBBD5", cex = 1.2)

abline(v = c(-logFCfilter, logFCfilter), lty = 2)
dev.off()

