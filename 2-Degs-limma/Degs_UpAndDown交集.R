# 加载包
library("VennDiagram")

# 设置工作目录（建议改为相对路径或通过参数传入）
setwd("D:/r/实验数据/2-Degs-limma/")

# 定义文件名列表
up_file_names <- c("GSE7158-upGenename.txt","GSE16134-upGenename.txt")
down_file_names <- c("GSE7158-downGenename.txt","GSE16134-downGenename.txt")
gene_lists_up <- list()
gene_lists_down <- list()
# 封装文件读取逻辑
read_gene_list <- function(file_name) {
  if (!file.exists(file_name)) {
    stop(paste("文件不存在:", file_name))
  }
  data <- tryCatch(
    read.table(file_name, header = FALSE, sep = "\t", check.names = FALSE),
    error = function(e) {
      stop(paste("读取文件失败:", file_name, "\n错误信息:", e$message))
    }
  )
  if (nrow(data) == 0) {
    warning(paste("文件为空:", file_name))
    return(character(0))
  }
  return(data[, 1])
}

# 批量读取基因列表
for (i in seq_along(up_file_names)) {
  gene_lists_up[[paste0("GSE", substr(up_file_names[i], 4, 9))]] <- read_gene_list(up_file_names[i])
}
# 批量读取基因列表
for (i in seq_along(down_file_names)) {
  gene_lists_down[[paste0("GSE", substr(down_file_names[i], 4, 9))]] <- read_gene_list(down_file_names[i])
}
# 绘制韦恩图（修复语法错误）
venn.diagram(
  x = gene_lists_up,
  filename = "VN-up.png",
  fill = rainbow(2),
  scaled = FALSE
)
# 绘制韦恩图（修复语法错误）
venn.diagram(
  x = gene_lists_down,
  filename = "VN-donw.png",
  fill = rainbow(2),
  scaled = FALSE
)
# 计算多个向量的交集（使用 Reduce 提高性能）
intersect_genes_up <- Reduce(intersect, gene_lists_up)
intersect_genes_down <- Reduce(intersect, gene_lists_down)
# 导出结果
if (length(intersect_genes_up) > 0) {
  write.table(
    file = "intersectGenes-up.txt",
    intersect_genes_up,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
} else {
  message("没有找到共同基因。")
}
# 导出结果
if (length(intersect_genes_down) > 0) {
  write.table(
    file = "intersectGenes-down.txt",
    intersect_genes_down,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
} else {
  message("没有找到共同基因。")
}

