# 加载包
library("VennDiagram")

# 设置工作目录（建议改为相对路径或通过参数传入）
setwd("G:/Rcode/实验数据/实验数据/2-Degs-limma")

# 定义文件名列表
file_names <- c("GSE10334-diffGenename.txt",  "GSE7158-diffGenename.txt")
gene_lists <- list()

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
for (i in seq_along(file_names)) {
  gene_lists[[paste0("GSE", substr(file_names[i], 4, 9))]] <- read_gene_list(file_names[i])
}

# 绘制韦恩图（修复语法错误）
venn.diagram(
  x = gene_lists,
  filename = "VN.png",
  fill = rainbow(2),
  scaled = FALSE
)

# 计算多个向量的交集（使用 Reduce 提高性能）
intersect_genes <- Reduce(intersect, gene_lists)

# 导出结果
if (length(intersect_genes) > 0) {
  write.table(
    file = "intersectGenes.txt",
    intersect_genes,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
} else {
  message("没有找到共同基因。")
}

