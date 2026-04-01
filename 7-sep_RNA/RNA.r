#############################################
# 🔥 Figure：Gingival tissues UMAP + MS4A6A
# 数据：GSE164241（牙龈单细胞）
#############################################
install.packages("Seurat")
install.packages("patchwork")
# ===== 1. 加载包 =====
library(Seurat)
library(ggplot2)
library(patchwork)

# ===== 2. 设置随机种子（保证UMAP稳定）=====
set.seed(123)
setwd("G:/Rcode/code/7-sep_RNA")

# ===== 3. 读取数据 =====
# 👉 修改为你的数据路径
data_dir <- "G:/Rcode/code/7-sep_RNA/GSE164241_RAW"

data <- Read10X(data.dir = data_dir)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = data,
  project = "Gingival_scRNA",
  min.cells = 3,
  min.features = 200
)

# ===== 4. 数据预处理 =====
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# ===== 5. 降维分析 =====
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:20
)

seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:20
)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.5
)

markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
pdf(paste0("MS4A6A_Umap.pdf"), width = 8, height = 6)
FeaturePlot(seurat_obj, features = c("MS4A6A"))
dev.off()
# ===== 6. 加载细胞注释 =====
# 👉 如果你有metadata文件（推荐）
# metadata <- read.csv("metadata.csv", row.names = 1)
# seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

# 👉 如果已有cell_type列
#Idents(seurat_obj) <- "cell_type"

# 如果没有注释，可临时用cluster代替：
 Idents(seurat_obj) <- "seurat_clusters"

# ===== 7. 右图：细胞类型UMAP =====
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
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

# ===== 8. 左图：MS4A6A表达 =====
p1 <- FeaturePlot(
  seurat_obj,
  features = "MS4A6A",
  reduction = "umap",
  pt.size = 0.5
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "red"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# ===== 9. 添加虚线框 + 标注 =====
# ⚠️ 需要根据你的UMAP实际调整坐标
p1 <- p1 +
  annotate("rect",
           xmin = -5, xmax = 1,
           ymin = -8, ymax = -2,
           linetype = "dashed",
           color = "black",
           fill = NA,
           size = 0.8) +
  annotate("text",
           x = -2,
           y = -1,
           label = "MS4A6A",
           size = 6)

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


