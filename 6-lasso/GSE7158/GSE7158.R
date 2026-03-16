# ==============================
# 0) 环境准备
# ==============================
rm(list = ls())

pkgs <- c("Boruta", "glmnet", "pROC", "data.table")
for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ==============================
# 1) 路径与输入文件（你按实际修改）
# ==============================
# 工作目录：放 shared_genes_list.txt + normalize表达矩阵 + Control/Treat.txt 的地方
setwd("D:/r/R_code/6-lasso/GSE7158")   # <- 改成你的训练集文件夹

# 训练集（PD）：表达矩阵 + 分组文件
train_expr_file <- "GSE7158-normalize.txt"
train_control_file <- "Control.txt"
train_treat_file   <- "Treat.txt"

# 外部验证集（PD）：通常放在另一个文件夹里
valid_dir <- "D:/r/R_code/2-Degs-limma/GSE56814"   # <- 改成你的验证集文件夹
valid_expr_file <- "GSE56814-normalize.txt"
valid_control_file <- "Control.txt"
valid_treat_file   <- "Treat.txt"

# shared基因列表：建议放在训练集目录下
shared_file <- "shared_genes_list.txt"

# 输出文件前缀
out_prefix <- "PD_ML_Boruta_LASSO"

# ==============================
# 2) 工具函数：读表达矩阵 + 对齐样本 + 构建 x/y
# ==============================
read_expr <- function(expr_file){
  df <- data.table::fread(expr_file, data.table = FALSE)
  # 兼容你的格式：第一列是ID
  if("ID" %in% colnames(df)){
    rownames(df) <- df$ID
    df$ID <- NULL
  } else {
    rownames(df) <- df[,1]
    df <- df[,-1]
  }
  m <- as.matrix(df)
  mode(m) <- "numeric"
  return(m)
}

read_group <- function(control_file, treat_file){
  ctrl <- read.table(control_file, header = FALSE, stringsAsFactors = FALSE)[,1]
  trt  <- read.table(treat_file,   header = FALSE, stringsAsFactors = FALSE)[,1]
  list(ctrl = ctrl, trt = trt)
}

build_xy <- function(expr_mat, ctrl_samples, trt_samples, genes_use){
  # 样本存在性检查
  stopifnot(all(ctrl_samples %in% colnames(expr_mat)))
  stopifnot(all(trt_samples  %in% colnames(expr_mat)))
  
  expr_sub <- expr_mat[, c(ctrl_samples, trt_samples), drop = FALSE]
  y <- factor(c(rep(0, length(ctrl_samples)), rep(1, length(trt_samples))))  # 0=Control, 1=Case
  
  g <- intersect(genes_use, rownames(expr_sub))
  if(length(g) < 3){
    stop(paste0("可用基因太少：", length(g), "。请检查 shared_genes 与表达矩阵行名是否一致（Symbol/探针ID）。"))
  }
  
  x <- t(expr_sub[g, , drop = FALSE])   # 样本×基因
  return(list(x = x, y = y, genes = g))
}

# ==============================
# 3) 读取 shared genes
# ==============================
shared_genes <- read.table(shared_file, header = FALSE, stringsAsFactors = FALSE)[,1]
shared_genes <- unique(shared_genes)
cat("Shared genes:", length(shared_genes), "\n")

# ==============================
# 4) 训练集：构建 x/y
# ==============================
train_expr <- read_expr(train_expr_file)
grp_train <- read_group(train_control_file, train_treat_file)

train_xy <- build_xy(train_expr, grp_train$ctrl, grp_train$trt, shared_genes)
x_train <- train_xy$x
y_train <- train_xy$y

cat("Train samples:", nrow(x_train), " | Train genes used:", ncol(x_train), "\n")

# ==============================
# 5) Boruta：特征筛选
# ==============================
set.seed(1)
bor <- Boruta(x = x_train, y = y_train, doTrace = 1)
bor_genes <- getSelectedAttributes(bor, withTentative = FALSE)

# 如果 Boruta 太严格筛到很少，可以把 tentative 也加上
if(length(bor_genes) < 5){
  bor_genes <- getSelectedAttributes(bor, withTentative = TRUE)
  cat("Boruta confirmed太少，已包含 tentative 特征。\n")
}

cat("Boruta genes:", length(bor_genes), "\n")
write.table(bor_genes, paste0(out_prefix, "_Boruta_genes.txt"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# ==============================
# 6) LASSO：压缩成签名基因 + 训练模型
# ==============================
x_bor <- as.matrix(x_train[, bor_genes, drop = FALSE])

set.seed(1)
fit_cv <- cv.glmnet(x_bor, y_train, alpha = 1, family = "binomial", nfolds = 10)

coef_min <- coef(fit_cv, s = "lambda.min")
lasso_genes <- rownames(coef_min)[which(as.numeric(coef_min) != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

cat("LASSO signature genes:", length(lasso_genes), "\n")
write.table(lasso_genes, paste0(out_prefix, "_LASSO_signature_genes.txt"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# ==============================
# 7) 训练集评估：AUC
# ==============================
pred_train <- as.numeric(predict(fit_cv, newx = x_bor, s = "lambda.min", type = "response"))
roc_train <- pROC::roc(y_train, pred_train, quiet = TRUE)
auc_train <- as.numeric(pROC::auc(roc_train))
cat("Train AUC =", auc_train, "\n")

# 保存ROC数据
train_roc_df <- data.frame(
  specificity = rev(roc_train$specificities),
  sensitivity = rev(roc_train$sensitivities)
)
write.csv(train_roc_df, paste0(out_prefix, "_Train_ROC_points.csv"), row.names = FALSE)
# ==============================
# 7.5) 绘制训练集 ROC 曲线
# ==============================
pdf(paste0(out_prefix, "_Train_ROC_Curve.pdf"), width = 6, height = 5)
plot(roc_train, 
     main = paste0("Training Set ROC Curve (AUC = ", round(auc_train, 3), ")"),
     col = "blue", 
     lwd = 2,
     print.auc = FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray")  # 添加对角线参考
dev.off()

png(paste0(out_prefix, "_Train_ROC_Curve.png"), width = 800, height = 600, res = 150)
plot(roc_train, 
     main = paste0("Training Set ROC Curve (AUC = ", round(auc_train, 3), ")"),
     col = "blue", 
     lwd = 2,
     print.auc = FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
# ==============================
# 8) 外部验证集：读入表达矩阵与分组，预测并算AUC
# ==============================
# 外部集可能在另一个目录
valid_expr <- NULL
grp_valid <- NULL

oldwd <- getwd()
setwd(valid_dir)
valid_expr <- read_expr(valid_expr_file)
grp_valid <- read_group(valid_control_file, valid_treat_file)
setwd(oldwd)

# 外部集使用同样的“训练基因空间”：Boruta genes（或 LASSO genes）
# 推荐用 Boruta genes 对齐，否则可能外部集少掉某些 LASSO gene导致无法预测
valid_xy <- build_xy(valid_expr, grp_valid$ctrl, grp_valid$trt, bor_genes)
x_valid <- valid_xy$x
y_valid <- valid_xy$y

# 对齐列顺序到训练时的 bor_genes
common_cols <- intersect(colnames(x_bor), colnames(x_valid))
x_valid2 <- as.matrix(x_valid[, common_cols, drop = FALSE])

# 同样把训练矩阵也裁剪到 common_cols，确保模型输入一致
x_train2 <- as.matrix(x_bor[, common_cols, drop = FALSE])

pred_valid <- as.numeric(predict(fit_cv, newx = x_valid2, s = "lambda.min", type = "response"))
roc_valid <- pROC::roc(y_valid, pred_valid, quiet = TRUE)
auc_valid <- as.numeric(pROC::auc(roc_valid))
cat("Valid AUC =", auc_valid, "\n")

valid_roc_df <- data.frame(
  specificity = rev(roc_valid$specificities),
  sensitivity = rev(roc_valid$sensitivities)
)
write.csv(valid_roc_df, paste0(out_prefix, "_Valid_ROC_points.csv"), row.names = FALSE)
# ==============================
# 8.5) 绘制验证集 ROC 曲线
# ==============================
pdf(paste0(out_prefix, "_Valid_ROC_Curve.pdf"), width = 6, height = 5)
plot(roc_valid, 
     main = paste0("Validation Set ROC Curve (AUC = ", round(auc_valid, 3), ")"),
     col = "red", 
     lwd = 2,
     print.auc = FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

png(paste0(out_prefix, "_Valid_ROC_Curve.png"), width = 800, height = 600, res = 150)
plot(roc_valid, 
     main = paste0("Validation Set ROC Curve (AUC = ", round(auc_valid, 3), ")"),
     col = "red", 
     lwd = 2,
     print.auc = FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
# ==============================
# 9) 保存汇总
# ==============================
summary_out <- data.frame(
  Train_AUC = auc_train,
  Valid_AUC = auc_valid,
  Shared_genes_used = ncol(x_train),
  Boruta_genes = length(bor_genes),
  LASSO_signature_genes = length(lasso_genes)
)
write.csv(summary_out, paste0(out_prefix, "_Summary.csv"), row.names = FALSE)

saveRDS(fit_cv, paste0(out_prefix, "_cvglmnet_model.rds"))
cat("全部完成 ✅ 结果已输出到当前工作目录。\n")

# ===== LASSO 两张标准图：系数路径 + CV曲线（带回执）=====
# 依赖：glmnet
# fit_cv <- cv.glmnet(x, y, alpha=1, family="binomial", nfolds=10)

# 1) 保存为 PDF（推荐）
# ===== 绘制带基因名的 LASSO 路径图 =====
pdf("LASSO_path_with_gene_names.pdf", width=7, height=5)

plot(fit_cv$glmnet.fit, xvar = "lambda", label = FALSE,
     main = "LASSO coefficient path",
     xlab = "log(Lambda)", ylab = "Coefficients")

# 取 lambda.min 下的系数（包含 Intercept）
coef_min <- coef(fit_cv, s = "lambda.min")

# 2) 计算“最右端”的 x 坐标（就用图的右边界）
usr <- par("usr")          # c(xmin, xmax, ymin, ymax)
x_right <- usr[2]          # 最右侧
x_text  <- x_right + 0.02*(usr[2]-usr[1])-0.4 # 稍微往右挪一点点（像你红框那样）
par(xpd=TRUE)              # 允许画到图框外一点点

# 3) 取每条曲线在最右端（最大lambda处）的系数
# beta: (p × nlambda)，最后一列就是最右端
beta_last <- as.matrix(fit_cv$glmnet.fit$beta)[, ncol(fit_cv$glmnet.fit$beta)]

# 4) 给每条曲线标基因名：只标“非零或绝对值较大”的，避免全部挤在0线上
# 你可以把阈值调小/调大
thr <- 1e-3
genes_to_label <- names(beta_last)[abs(beta_last) > thr]

for(g in genes_to_label){
  y_pos <- beta_last[g]
  text(x = x_text, y = y_pos, labels = g, pos = 4, cex = 0.6, col = "black")
}
# 画 lambda.min / lambda.1se 竖线（建议加上，更像论文标准图）
# abline(v = log(fit_cv$lambda.min), lty = 2, col = "red")
# abline(v = log(fit_cv$lambda.1se), lty = 3, col = "blue")

legend("topright",
       legend = c(paste0("lambda.min=", signif(fit_cv$lambda.min, 3)),
                  paste0("lambda.1se=", signif(fit_cv$lambda.1se, 3))),
       lty = c(2, 3), col = c("red", "blue"),
       bty = "n", cex = 0.8)

dev.off()
# 绘制 LASSO 模型系数路径图和 CV 曲线
pdf("LASSO_path_and_CV2.pdf", width=6, height=8)
# (B) 交叉验证曲线
plot(fit_cv, main="Cross-validation curve",
     xlab="log(Lambda)", ylab="Binomial deviance")
abline(v=log(fit_cv$lambda.min), lty=2)
abline(v=log(fit_cv$lambda.1se), lty=3)

dev.off()

# 2) 同时保存 PNG（可选）
png("LASSO_path_and_CV.png", width=900, height=1200, res=150)
par(mfrow=c(2,1), mar=c(4,4,3,1))

plot(fit_cv$glmnet.fit, xvar = "lambda", label = FALSE,
     main = "LASSO coefficient path",
     xlab = "log(Lambda)", ylab = "Coefficients")

# 取 lambda.min 下的系数（包含 Intercept）
coef_min <- coef(fit_cv, s = "lambda.min")

# 右侧x坐标（log(lambda) 最大的位置）
x_pos <- max(log(fit_cv$glmnet.fit$lambda))
print(x_pos)

# 逐个基因标注（按基因名取系数，避免错位）
for(g in lasso_genes){
  if(g %in% rownames(coef_min)){
    y_pos <- as.numeric(coef_min[g, 1])
    # 只标非零系数的（可选，避免一堆0挤在0线上）
    if(abs(y_pos) > 1e-6){
      text(x = x_pos, y = y_pos, labels = g,
           pos = 4, cex = 0.6, col = "black", xpd = TRUE)
    }
  }
}

# 画 lambda.min / lambda.1se 竖线（建议加上，更像论文标准图）
# abline(v = log(fit_cv$lambda.min), lty = 2, col = "red")
# abline(v = log(fit_cv$lambda.1se), lty = 3, col = "blue")

legend("topright",
       legend = c(paste0("lambda.min=", signif(fit_cv$lambda.min, 3)),
                  paste0("lambda.1se=", signif(fit_cv$lambda.1se, 3))),
       lty = c(2, 3), col = c("red", "blue"),
       bty = "n", cex = 0.8)

plot(fit_cv, main="Cross-validation curve",
     xlab="log(Lambda)", ylab="Binomial deviance")
abline(v=log(fit_cv$lambda.min), lty=2)
abline(v=log(fit_cv$lambda.1se), lty=3)
dev.off()

# ==============================
# 9.5) 绘制训练集与验证集对比 ROC 图
# ==============================
pdf(paste0(out_prefix, "_Compare_ROC_Curve.pdf"), width = 6, height = 5)
plot(roc_train, 
     col = "blue", 
     lwd = 2,
     main = "ROC Curves: Training vs Validation",
     print.auc = FALSE)
plot(roc_valid, 
     col = "red", 
     lwd = 2, 
     add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", 
       legend = c(paste0("Train AUC=", round(auc_train, 3)), 
                  paste0("Valid AUC=", round(auc_valid, 3))),
       col = c("blue", "red"), 
       lwd = 2,
       bty = "n")
dev.off()

png(paste0(out_prefix, "_Compare_ROC_Curve.png"), width = 800, height = 600, res = 150)
plot(roc_train, 
     col = "blue", 
     lwd = 2,
     main = "ROC Curves: Training vs Validation",
     print.auc = FALSE)
plot(roc_valid, 
     col = "red", 
     lwd = 2, 
     add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", 
       legend = c(paste0("Train AUC=", round(auc_train, 3)), 
                  paste0("Valid AUC=", round(auc_valid, 3))),
       col = c("blue", "red"), 
       lwd = 2,
       bty = "n")
dev.off()

