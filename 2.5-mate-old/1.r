# ===== 0) 安装/加载包 =====
if (!require(metafor)) install.packages("metafor")
library(metafor)
setwd("D:/r/实验数据/2.5-mate")
# ===== 1) 读取 limma 全基因结果的函数（兼容你的 id 列）=====
read_limma_allgene <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  # 兼容：如果第一列不是id但文件里有行名型基因名
  if (!("id" %in% colnames(df))) {
    # 常见情况：第一列是基因名，但列名不是id
    colnames(df)[1] <- "id"
  }
  
  # 兼容：如果文件第一行是你之前写入的“id\tlogFC...”那种重复表头
  # 表现为：df$id 里出现 "logFC"/"AveExpr"/"t" 等关键词
  bad <- df$id %in% c("logFC","AveExpr","t","P.Value","adj.P.Val","B","id")
  if (any(bad)) df <- df[!bad, ]
  
  # 确保数值列是数值
  num_cols <- intersect(c("logFC","t","P.Value","adj.P.Val"), colnames(df))
  for (cc in num_cols) df[[cc]] <- as.numeric(df[[cc]])
  
  # 去重：同一基因多行时保留最小P.Value那行
  df <- df[!is.na(df$id) & df$id != "", ]
  df <- df[order(df$P.Value), ]
  df <- df[!duplicated(df$id), ]
  rownames(df) <- df$id
  df
}

# ===== 2) 读入你的四个结果表 =====
pd1 <- read_limma_allgene("GSE10334-all.gene.txt")
pd2 <- read_limma_allgene("GSE16134-all.gene.txt")

op1 <- read_limma_allgene("GSE56814-all.gene.txt")
op2 <- read_limma_allgene("GSE56815-all.gene.txt")
# ===== 3) PD meta：合并 logFC（需要 SE，用 SE = abs(logFC/t) 估计）=====
pd_common <- intersect(rownames(pd1), rownames(pd2))
pd1c <- pd1[pd_common, ]
pd2c <- pd2[pd_common, ]

pd1c$SE <- abs(pd1c$logFC / pd1c$t)
pd2c$SE <- abs(pd2c$logFC / pd2c$t)

pd_long <- rbind(
  data.frame(id = rownames(pd1c), yi = pd1c$logFC, sei = pd1c$SE, study = "PD1"),
  data.frame(id = rownames(pd2c), yi = pd2c$logFC, sei = pd2c$SE, study = "PD2")
)
pd_long <- pd_long[is.finite(pd_long$sei) & pd_long$sei > 0, ]

pd_meta <- do.call(rbind, lapply(split(pd_long, pd_long$id), function(d) {
  if (nrow(d) < 2) return(NULL)
  fit <- rma.uni(yi = d$yi, sei = d$sei, method = "FE")
  data.frame(
    id = d$id[1],
    PD_meta_logFC = as.numeric(fit$b),
    PD_meta_p = as.numeric(fit$pval),
    PD_I2 = as.numeric(fit$I2)
  )
}))

pd_meta$PD_meta_FDR <- p.adjust(pd_meta$PD_meta_p, method = "BH")


write.table(pd_meta, "PD_meta_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# ===== 4) OP meta：Fisher 合并 p 值 + 方向一致性 =====
op_common <- intersect(rownames(op1), rownames(op2))
op1c <- op1[op_common, ]
op2c <- op2[op_common, ]

# Fisher: X = -2*sum(log(p)) ~ chisq(df=2k)
op_meta_p <- pchisq(-2 * (log(op1c$P.Value) + log(op2c$P.Value)), df = 4, lower.tail = FALSE)

# 方向一致性：两队列 logFC 同号
dir_ok <- sign(op1c$logFC) == sign(op2c$logFC) & sign(op1c$logFC) != 0
op_meta <- data.frame(
  id = op_common,
  OP_meta_p = op_meta_p,
  OP_dir_consistent = dir_ok,
  OP_mean_logFC = (op1c$logFC + op2c$logFC) / 2
)
op_meta$OP_meta_FDR <- p.adjust(op_meta$OP_meta_p, method = "BH")

write.table(op_meta, "OP_meta_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# ===== 5) shared genes =====
pd_sig <- subset(pd_meta, PD_meta_FDR < 0.1)
op_sig <- subset(op_meta, OP_meta_FDR < 0.1 & OP_dir_consistent)

shared_ids <- intersect(pd_sig$id, op_sig$id)
set.seed(123456)
# 方向一致：PD_meta_logFC 与 OP_mean_logFC 同号
tmp <- merge(pd_sig[, c("id","PD_meta_logFC","PD_meta_FDR")],
             op_sig[, c("id","OP_mean_logFC","OP_meta_FDR")],
             by = "id", all = FALSE)

shared_final <- subset(tmp, sign(PD_meta_logFC) == sign(OP_mean_logFC) & sign(PD_meta_logFC) != 0)

write.table(shared_final, "shared_genes_meta.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(shared_final$id, "shared_genes_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat("shared genes 数量 =", nrow(shared_final), "\n")
# ===== 6) ML：以PD训练集为例（需要表达矩阵）=====
if (!require(Boruta)) install.packages("Boruta")
if (!require(glmnet)) install.packages("glmnet")
if (!require(pROC)) install.packages("pROC")
library(Boruta); library(glmnet); library(pROC)

shared_genes <- read.table("shared_genes_list.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

# 读取训练集表达矩阵（你自己替换文件名）
train_expr <- read.table("GSE56815-normalize.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(train_expr) <- train_expr$ID
train_expr$ID <- NULL
train_expr <- as.matrix(train_expr); mode(train_expr) <- "numeric"

# 训练集分组（由Control.txt Treat.txt来）
Control <- read.table("Control.txt", header=FALSE, stringsAsFactors=FALSE)[,1]
Treat <- read.table("Treat.txt", header=FALSE, stringsAsFactors=FALSE)[,1]
train_expr <- train_expr[, c(Control, Treat), drop=FALSE]
y <- factor(c(rep(0, length(Control)), rep(1, length(Treat))))  # 0=Control, 1=Case

# 取 shared genes（注意：基因名要能对上）
genes_use <- intersect(shared_genes, rownames(train_expr))
x <- t(train_expr[genes_use, , drop=FALSE])  # 样本×基因

# Boruta
set.seed(1)
bor <- Boruta(x = x, y = y, doTrace = 1)
bor_genes <- getSelectedAttributes(bor, withTentative = FALSE)

# LASSO
x2 <- as.matrix(x[, bor_genes, drop=FALSE])
fit <- cv.glmnet(x2, y, alpha = 1, family = "binomial")
coef_min <- coef(fit, s = "lambda.min")
lasso_genes <- rownames(coef_min)[which(as.numeric(coef_min) != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

cat("Boruta genes:", length(bor_genes), " | LASSO genes:", length(lasso_genes), "\n")
write.table(lasso_genes, "ML_signature_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# 训练集预测与AUC
pred <- as.numeric(predict(fit, newx = x2, s = "lambda.min", type = "response"))
auc_train <- auc(roc(y, pred))
cat("Train AUC =", as.numeric(auc_train), "\n")

