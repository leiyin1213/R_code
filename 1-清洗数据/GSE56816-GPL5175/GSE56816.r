# ===== 0) packages =====
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(GEOquery)
library(stringr)
library(dplyr)

GSE <- "GSE56816"  #数据集名称
GPL<- "GPL5175-3188.txt" #芯片名称
geneName<-"gene_assignment" #芯片对应的基因行名
base_path <- "D:/r/实验数据/1-清洗数据/GSE56816-GPL5175"  #数据存放路径
#full_path <- file.path(base_path, GSE)
setwd(base_path)

# ===== 1) download GEO =====

gset <- getGEO(GSE,destdir = base_path,AnnotGPL= F, getGPL= F)

dat = exprs(gset[[1]])#使用GPL5175芯片版本

ex <- dat
qx <- as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm = T))

logC <-(qx[5]>100 || 
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2]<1 && qx[4]>1&& qx[4]<2 )
)
if(logC){ 
  ex[which(ex <= 0)] <- NA
  dat <- log2(ex)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
#获取临床信息
pd <- pData(gset[[1]])

write.csv(pd,paste0("Clinical_", GSE,".csv"),row.names = TRUE)
#读取txt注释文件
gpl= read.table(GPL,
                header = TRUE ,
                fill= T ,
                sep="\t",
                comment.char = "#",
                stringsAsFactors=FALSE,
                quote = "")
#查看
#View(gpl)
colnames(gpl)

#提取探针Id及基因symbol
ids= gpl[,c("ID",geneName)]

#修改列名
colnames(ids)=c('probe_id',"Symbol")

#获取基因symbol
library(stringr)
ids$Symbol= trimws(str_split(ids$`Symbol`,"//",simplify= T)[,2]) #根据符号进行分割(需要按照芯片数据修改)
#去除没有注释的探针
ids= ids[ids$Symbol !='',]
ids= ids[ids$Symbol !='---',]
ids= ids[ids$Symbol !='---',]
#平台文件的ID和矩阵种的ID匹配，%in%用于判断是否匹配
ids= ids[as.character(ids$probe_id) %in%as.character(rownames(dat)),]

#匹配表达数据
dat= dat[as.character(ids$probe_id),]

#查看探针名于注释文件名是否一致
table(as.character(rownames(dat))==as.character(ids$probe_id))

#合并
dat <- cbind(ids,dat)

#删除重复基因

#dat <-aggregate( .~Symbol,data=dat,max)
dat <- dat %>%
  as.data.frame() %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")#取同基因平均值
# 1. 强制转换为基础 data.frame
dat <- as.data.frame(dat)

# 2. 提取第一列作为向量 (使用 [[1]] 而不是 [,1])
row_names <- dat[[1]]

# 3. 设置行名
rownames(dat) <- row_names

# 4. 删除第一列 (因为已经变为行名)
dat <- dat[,-1]
dat <- dat[,-1]
write.table(data.frame(ID=rownames(dat),dat),file=paste0(GSE,".txt"),sep = '\t',quote=F,row.names = F)

