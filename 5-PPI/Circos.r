install.packages(c("circlize", "dplyr"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("biomaRt"))


hub <- c("MAP1LC3B",
         "TUBA8",
         "RAB2A",
         "DAPK2",
         "DAP3",
         "CDYL",
         "SFSWAP",
         "UTP6",
         "PTMA",
         "PFDN5")   # 这里换成你自己的hub基因

library(biomaRt)
library(dplyr)

mart <- useMart(host = "https://useast.ensembl.org/", 
                biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")

loc <- getBM(
  attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"),
  filters = "hgnc_symbol",
  values = hub,
  mart = mart
)

loc <- loc %>%
  filter(chromosome_name %in% c(1:22,"X","Y")) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  select(gene=hgnc_symbol, chr, start=start_position, end=end_position)


library(circlize)

circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.01, 0.01))

# 初始化染色体扇区：需要每条染色体长度
# 用 biomaRt 再抓一次 chr length（简单起见也可以用内置长度，但我这里给你可复制写法）
chr_len <- getBM(
  attributes = c("listAttributes", "length"),
  filters = "listAttributes",
  values = c(as.character(1:22),"X","Y"),
  mart = mart
) %>% mutate(chr=paste0("chr", listAttributes)) %>% select(chr, length)

circos.initialize(factors = chr_len$chr, xlim = cbind(rep(0,nrow(chr_len)), chr_len$length))

# 画染色体轨道
circos.trackPlotRegion(ylim = c(0,1), bg.border = NA, track.height = 0.06,
                       panel.fun = function(x, y) {
                         chr = CELL_META$sector.index
                         circos.text(CELL_META$xcenter, 1.2, gsub("chr","",chr), cex=0.5, facing="clockwise", niceFacing=TRUE)
                       }
)

# 把基因位置画进去（点 + 基因名）
circos.trackPlotRegion(ylim=c(0,1), bg.border = NA, track.height = 0.10,
                       panel.fun=function(x,y){ }
)

apply(loc, 1, function(r){
  circos.points(r["start"], 0.5, sector.index = r["chr"], pch=16, cex=0.6)
  circos.text(as.numeric(r["start"]), 0.8, r["gene"], sector.index=r["chr"],
              cex=0.45, facing="clockwise", niceFacing=TRUE)
})