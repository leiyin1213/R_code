BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scrapper")

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(stringr)
library(SingleR)
library(celldex)
library(scrapper)
library(dplyr)
library(ggsci)  # 添加此行
library(scales)
# =========================
# 🔧 1. 设置路径
# =========================
setwd("G:/Rcode/code/7-sep_RNA")
base_dir <- "G:/Rcode/code/7-sep_RNA/GSE164241_RAW"
hpca.se <- HumanPrimaryCellAtlasData()
ref2 <- BlueprintEncodeData()
# =========================
# 🔍 2. 找所有 matrix 文件（每个代表一个样本）
# =========================
all_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)

mtx_files <- all_files[grepl("matrix.mtx", all_files)]

if(length(mtx_files) == 0){
  stop("❌ 没找到 matrix.mtx 文件")
}

cat("✅ 找到", length(mtx_files), "个样本\n")

# =========================
# 📦 3. 逐个读取样本
# =========================

seurat_list <- list()

for(i in seq_along(mtx_files)){
  
  mtx_path <- mtx_files[i]
  
  cat("\n🚀 处理:", mtx_path, "\n")
  
  # 🔥 提取 GSM ID（关键！！）
  gsm_id <- str_extract(basename(mtx_path), "GSM[0-9]+")
  
  cat("匹配样本:", gsm_id, "\n")
  
  # 🔥 精准匹配对应文件
  feature_file <- all_files[grepl(gsm_id, all_files) & grepl("features|genes", all_files)]
  barcode_file <- all_files[grepl(gsm_id, all_files) & grepl("barcodes", all_files)]
  
  if(length(feature_file) == 0 | length(barcode_file) == 0){
    stop(paste("❌ 找不到匹配文件:", gsm_id))
  }
  
  # 读取
  mat <- readMM(mtx_path)
  features <- read.delim(feature_file[1], header = FALSE)
  barcodes <- read.delim(barcode_file[1], header = FALSE)
  
  # 检查（非常重要）
  cat("matrix:", dim(mat), "\n")
  cat("features:", nrow(features), "\n")
  cat("barcodes:", nrow(barcodes), "\n")
  
  # ❗ 严格检查（不再自动截断）
  if(nrow(features) != nrow(mat)){
    stop(paste("❌ 基因数不匹配:", gsm_id))
  }
  
  if(nrow(barcodes) != ncol(mat)){
    stop(paste("❌ 细胞数不匹配:", gsm_id))
  }
  
  # 设置名称
  gene_names <- if(ncol(features) >= 2) features[,2] else features[,1]
  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- paste0(gsm_id, "_", barcodes[,1])
  
  # 创建对象
  obj <- CreateSeuratObject(mat, project = gsm_id)
  obj$sample <- gsm_id
  
  seurat_list[[gsm_id]] <- obj
  
  cat("✅ 完成:", gsm_id, "\n")
}

# =========================
# 🔗 合并
# =========================
combined <- merge(seurat_list[[1]], y = seurat_list[-1])
combined <- JoinLayers(combined)
cat("🎉 合并完成！细胞数:", ncol(combined), "\n")

# =========================
# 🔍 5. QC
# =========================
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

VlnPlot(combined, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample",
        ncol = 3)

# 过滤
combined <- subset(combined,
                   subset = nFeature_RNA > 300 &
                            nFeature_RNA < 6000 &
                            percent.mt < 10)

cat("✅ QC完成\n")

# =========================
# ⚖️ 6. 标准流程
# =========================
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)

combined <- ScaleData(combined)
combined <- RunPCA(combined)

ElbowPlot(combined)

# =========================
# 🌌 7. UMAP + 聚类
# =========================
combined <- RunUMAP(combined, dims = 1:20)

combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# 看聚类
DimPlot(combined, label = TRUE)

# 看样本来源（非常重要！）
DimPlot(combined, group.by = "sample")

cat("🎉 UMAP完成！\n")

# =========================
# 🏷️ 8. Marker
# =========================
markers <- FindAllMarkers(combined,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

head(markers)
saveRDS(markers, file = "GSE164241_markers.rds")
# =========================
# 💾 保存
# =========================
saveRDS(combined, file = "GSE164241_combined.rds")

markers <- readRDS("GSE164241_markers.rds")
combined<-readRDS("GSE164241_combined.rds")
# =========================
# 💾 自动注释
# =========================

###把rna的转录表达数据提取
testdata <- GetAssayData(
  combined,
  assay = "RNA",
  layer  = "data"
)
clusters <- combined$seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref2, 
                    labels = ref2$label.main, 
                    clusters = clusters,assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
##添加到metadata当中
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
combined@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  combined@meta.data[which(combined@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

DimPlot(combined, reduction = "umap",label = T)
DimPlot(combined, reduction = "umap", group.by = "celltype",label = F)
#更改Active Idents
#重新划分细胞亚群
cluster_map <- combined@meta.data %>%
  group_by(seurat_clusters, celltype) %>%
  summarise(n = n()) %>%
  slice_max(n, n = 1)

new.cluster.ids <- cluster_map$celltype

names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap",label = F)
save(combined,file="OS.relabel.Rdata")
load("OS.relabel.Rdata")



# ===== 7. 右图：细胞类型UMAP =====
FeaturePlot(combined, features = c("MS4A6A"))
VlnPlot(combined,features = c("MS4A6A"),pt.size = 0)



p2 <- DimPlot(
  combined,
  reduction = "umap",
 # group.by = "seurat_clusters",
  pt.size = 0.5
) +
  theme_classic() +
  ggtitle("gingival tissues") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

p2
# ===== 8. 左图：MS4A6A表达 =====
p1 <- FeaturePlot(
  combined,
  reduction = "umap",
  features = c("MS4A6A"),
  pt.size = 0.5
) 
+
  scale_color_gradient(
    low = "lightgrey",
    high = "red"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )
# p1
# 获取细胞类型信息
cell_types <- unique(combined$celltype)
tmp_length <- length(cell_types)

# cell_type_cols <- scale_color_hue(tmp_length)
# 使用 RColorBrewer 或 rainbow 生成颜色
# 使用内置的颜色函数替代 brewer.pal
if(tmp_length <= 12) {
  # 使用 rainbow 或其他内置调色板
  cell_type_cols <- rainbow(max(tmp_length, 3))
  # 或者使用 gray.colors
  # cell_type_cols <- gray.colors(max(tmp_length, 3))
} else {
  # 如果细胞类型超过12个，使用彩虹色
  cell_type_cols <- rainbow(tmp_length)
}

# 找到特定细胞类型和细胞ID
epithelial_cells <- WhichCells(combined, idents = "Melanocytes")  # 根据实际的细胞类型ID调整

# 绘制高亮图
DimPlot(combined, 
        cells.highlight = list(Epithelial = epithelial_cells),
        label = TRUE,
        reduction = "umap",
        cols.highlight = "red",  # 高亮颜色
        cols = cell_type_cols)   # 背景色
# ===== 9. 添加虚线框 + 标注 =====
# ⚠️ 需要根据你的UMAP实际调整坐标
p1 <- p1 +
  annotate("rect",
           xmin = -2, xmax = 5,
           ymin = -8, ymax = -1,
           linetype = "dashed",
           color = "black",
           fill = NA,
           size = 0.8) +
  annotate("text",
           x = 2,
           y = 0,
           label = "MS4A6A",
           size = 6)
p1
# ===== 10. 拼图 =====
p_final <- p1 + p2

# 显示
print(p_final)

# ===== 11. 保存图片 =====
ggsave(
  filename = "Figure_MS4A6A_gingival_umap.png",
  plot = p_final,
  width = 12,
  height = 12,
  dpi = 300
)
cat("✅ 全流程完成！\n")

