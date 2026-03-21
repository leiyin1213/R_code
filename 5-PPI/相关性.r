# ===============================
# 0. 自动安装缺失包
# ===============================
packages <- c("corrplot")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ===============================
# 1. 读取表达矩阵
# ===============================

setwd("D:/r/R_code/5-PPI")


expr <- read.table("GSE7158.txt",
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

# 设置行名（基因）
rownames(expr) <- expr[,1]

# 删除前两列（ID + symbol）
expr <- expr[,-c(1,2)]

# 转为数值矩阵
expr <- as.matrix(expr)
mode(expr) <- "numeric"


# ===============================
# 2. 选择基因（和你图一致）
# ===============================
genes_use <-read.table("Top4_intersection_MCCsorted.txt",sep="\t",header=T,check.names=F) [,1]

# 只保留存在的基因
genes_use <- intersect(genes_use, rownames(expr))

expr_sub <- expr[genes_use, ]


# ===============================
# 3. 选一部分样本（先测试用）
# ===============================
expr_PD <- expr_sub[1:10, 1:10]


# ===============================
# 4. 计算相关性
# ===============================
cor_mat <- cor(t(expr_sub),
               method = "pearson",
               use = "pairwise.complete.obs")


# ===============================
# 5. 绘图
# ===============================
corrplot(cor_mat,
         method = "pie",
         type = "lower",
         
         # 颜色（偏蓝风格）
         col = colorRampPalette(c("red","white","#053061"))(200),
         
         # 圆更大
         tl.cex = 1,
         tl.col = "black",
         tl.srt = 90,   # 👈 竖直！
         
         # 数值
         addCoef.col = "black",
         number.cex = 0.6,
         
         # 色条
         cl.pos = "b",
         cl.lim = c(-1,1),
         
         # 去边框感
         bg = "white",
         
         # 标题
         title = "PD",
         mar = c(0,0,2,0)
)

# ===============================
# 6. 保存为PDF（可选）
# ===============================
pdf("PD_correlation_plot.pdf", width=6, height=6)

corrplot(cor_mat,
         method = "circle",
         type = "lower",
         tl.col = "black",
         tl.cex = 0.8,
         col = colorRampPalette(c("#8B0000","white","#00008B"))(200),
         addCoef.col = "black",
         number.cex = 0.7)

dev.off()

