install.packages(c("circlize", "dplyr"))
library(circlize)
library(dplyr)

# =========================
# 1. 手动输入你的基因坐标（完全离线）
# 这里先用你当前这4个基因
# =========================
loc <- data.frame(
  chr   = c("chr11", "chr1", "chr17", "chr12"),
  start = c(612876, 94884755, 81510341, 56729177),
  end   = c(614726, 94885509, 81523876, 56746337),
  gene  = c("IRF7", "ISG15", "IFI35", "STAT2"),
  stringsAsFactors = FALSE
)

# =========================
# 2. 人类染色体长度（hg38，完全离线）
# =========================
chr_len <- data.frame(
  chr = c(paste0("chr", 1:22), "chrX", "chrY"),
  length = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
    90338345, 83257441, 80373285, 58617616, 64444167,
    46709983, 50818468, 156040895, 57227415
  ),
  stringsAsFactors = FALSE
)

# =========================
# 3. 只给有基因的染色体加外圈颜色
# 你可以自己换颜色
# =========================
chr_col <- c(
  "chr1"  = "red",
  "chr11" = "forestgreen",
  "chr12" = "blue",
  "chr17" = "purple"
)

# =========================
# 4. 开始画图
# =========================
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 1,
  track.margin = c(0.002, 0.002),
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

circos.initialize(
  factors = chr_len$chr,
  xlim = cbind(rep(0, nrow(chr_len)), chr_len$length)
)

# =========================
# 5. 外圈彩色高亮线
# =========================
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.02,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    
    if (chr %in% names(chr_col)) {
      circos.segments(
        x0 = xlim[1], y0 = 0.5,
        x1 = xlim[2], y1 = 0.5,
        col = chr_col[chr],
        lwd = 3
      )
    }
  }
)

# =========================
# 6. 染色体主环（灰白相间，模拟参考图风格）
# 不是严格cytoband，但视觉上接近
# =========================
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    chr_length <- xlim[2]
    
    # 整条染色体底色
    circos.rect(
      xleft = xlim[1], ybottom = 0,
      xright = xlim[2], ytop = 1,
      col = "white", border = "black", lwd = 0.8
    )
    
    # 画一些交替的小条，模拟黑白核型带
    n_band <- 18
    band_breaks <- seq(0, chr_length, length.out = n_band + 1)
    
    for (i in seq_len(n_band)) {
      if (i %% 2 == 0) {
        circos.rect(
          xleft = band_breaks[i],
          ybottom = 0,
          xright = band_breaks[i + 1],
          ytop = 1,
          col = "grey75",
          border = NA
        )
      }
    }
  }
)

# =========================
# 7. 染色体编号
# =========================
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.10,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    circos.text(
      x = CELL_META$xcenter,
      y = 0.5,
      labels = gsub("chr", "", chr),
      cex = 0.7,
      facing = "clockwise",
      niceFacing = TRUE
    )
  }
)

# =========================
# 8. 基因引线和标签
# =========================
for (i in 1:nrow(loc)) {
  ch <- loc$chr[i]
  x  <- loc$start[i]
  gene_name <- loc$gene[i]
  
  # 从染色体环向内画短引线
  circos.segments(
    sector.index = ch,
    x0 = x, y0 = 0.78,
    x1 = x, y1 = 0.45,
    col = "grey40",
    lwd = 0.8
  )
  
  # 基因名字
  circos.text(
    sector.index = ch,
    x = x,
    y = 0.38,
    labels = gene_name,
    cex = 0.45,
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5)
  )
}