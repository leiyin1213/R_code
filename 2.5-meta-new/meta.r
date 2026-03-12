# ===== 0) 安装/加载包 =====
if (!require(metafor)) install.packages("metafor")
library(metafor)
setwd("D:/r/实验数据/2.5-meta-new")

# ===== 1) 读取 limma 全基因结果（要求包含 SE）=====
read_limma_allgene2 <- function(file){
  df <- read.table(file, header=TRUE, sep="\t", check.names=FALSE,
                   stringsAsFactors=FALSE, quote="", comment.char="")
  # 基因名在行名里（你现在文件就是这样）
  if (is.null(rownames(df)) || all(rownames(df) == as.character(1:nrow(df)))) {
    stop(paste0("No gene names in rownames for file: ", file))
  }
  # 强制数值列
  for(cc in intersect(c("logFC","AveExpr","t","P.Value","adj.P.Val","B","SE"), colnames(df))){
    df[[cc]] <- as.numeric(df[[cc]])
  }
  df$P.Value <- pmax(df$P.Value, 1e-300)
  # 去重（万一有重复基因名）
  df <- df[order(df$P.Value), ]
  df <- df[!duplicated(rownames(df)), ]
  df
}

# ===== 2) 读入四个结果表 =====
pd1 <- read_limma_allgene2("GSE10334-all.gene.txt")
pd2 <- read_limma_allgene2("GSE16134-all.gene.txt")
op1 <- read_limma_allgene2("GSE56814-all.gene.txt")
op2 <- read_limma_allgene2("GSE56815-all.gene.txt")

length(intersect(rownames(pd1), rownames(pd2)))  # 现在一定 > 0

# ===== 通用：对两队列做效应量 meta（每个基因一个 rma.uni）=====
run_meta_two_cohorts <- function(dfA, dfB, prefix="PD", method="FE") {
  common <- intersect(rownames(dfA), rownames(dfB))
  if(length(common) == 0) stop(paste0(prefix, ": no common gene IDs. Check rownames."))
  
  a <- dfA[common, , drop=FALSE]
  b <- dfB[common, , drop=FALSE]
  
  long <- rbind(
    data.frame(id = common, yi = a$logFC, sei = a$SE, study = "A"),
    data.frame(id = common, yi = b$logFC, sei = b$SE, study = "B")
  )
  
  long <- long[is.finite(long$sei) & long$sei > 0 & is.finite(long$yi), ]
  
  meta <- do.call(rbind, lapply(split(long, long$id), function(d) {
    fit <- rma.uni(yi = d$yi, sei = d$sei, method = method)
    data.frame(
      id = d$id[1],
      meta_logFC = as.numeric(fit$b),
      meta_se = as.numeric(fit$se),
      meta_p = as.numeric(fit$pval),
      I2 = if (!is.null(fit$I2)) as.numeric(fit$I2) else NA_real_,
      tau2 = if (!is.null(fit$tau2)) as.numeric(fit$tau2) else NA_real_
    )
  }))
  
  meta$meta_FDR <- p.adjust(meta$meta_p, method="BH")
  
  # 加前缀
  names(meta) <- sub("^meta_", paste0(prefix, "_meta_"), names(meta))
  names(meta) <- sub("^I2$", paste0(prefix, "_I2"), names(meta))
  names(meta) <- sub("^tau2$", paste0(prefix, "_tau2"), names(meta))
  meta
}

# ===== 3) PD meta：推荐先 FE（同平台），并输出I2供判断 =====
pd_meta <- run_meta_two_cohorts(pd1, pd2, prefix="PD", method="FE")      # 同平台可先FE
op_meta_FE <- run_meta_two_cohorts(op1, op2, prefix="OP", method="FE")
op_meta_RE <- run_meta_two_cohorts(op1, op2, prefix="OP", method="REML")  # 作为敏感性分析


# （可选）PD 也跑一份随机效应做敏感性分析
pd_meta_RE <- run_meta_two_cohorts(pd1, pd2, prefix="PD", method="REML")
write.table(pd_meta_RE, "PD_meta_REML.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ===== 5) shared genes：用 meta_FDR + 效应量阈值 + 方向一致 =====
logFC_thr <- log2(1.1)   # 0.263
fdr_thr <- 0.1           # 探索用0.1，成稿建议0.05

pd_sig <- subset(pd_meta, PD_meta_FDR < fdr_thr & abs(PD_meta_logFC) > logFC_thr)
op_sig <- subset(op_meta_FE, OP_meta_FDR < fdr_thr & abs(OP_meta_logFC) >logFC_thr )

common_op <- intersect(rownames(op1), rownames(op2))
dir_ok_op <- sign(op1[common_op,"logFC"]) == sign(op2[common_op,"logFC"]) & sign(op1[common_op,"logFC"])!=0
dir_ok_df <- data.frame(id = common_op, OP_dir_ok = dir_ok_op)

op_meta_FE2 <- merge(op_meta_FE, dir_ok_df, by="id", all.x=TRUE)
op_sig <- subset(op_meta_FE, OP_meta_FDR < fdr_thr & abs(OP_meta_logFC) >logFC_thr & OP_dir_ok)

shared <- merge(pd_sig[,c("id","PD_meta_logFC","PD_meta_FDR")],
                op_sig[,c("id","OP_meta_logFC","OP_meta_FDR")],
                by="id")

shared_final <- subset(shared, sign(PD_meta_logFC) == sign(OP_meta_logFC) & sign(PD_meta_logFC) != 0)

write.table(shared_final, "shared_genes_meta.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(shared_final$id, "shared_genes_list.txt", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

cat("shared genes 数量 =", nrow(shared_final), "\n")

