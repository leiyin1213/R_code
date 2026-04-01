library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)
# =========================
# 🔧 1. 设置路径
# =========================
base_dir <- "G:/Rcode/code/7-sep_RNA/GSE164241_RAW"

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

# =========================
# 💾 保存
# =========================
saveRDS(combined, file = "GSE164241_combined.rds")

cat("✅ 全流程完成！\n")