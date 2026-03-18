# 设置工作目录
setwd("/Volumes/Seagate Basic/精品代码/29_Lasso回归筛选特征基因4张图")
# 加载所需包
library(glmnet)  # 用于Lasso回归分析
# 设置输入文件
input_data <- "training_data.csv"  
# 读取输入数据
expression_data <- read.csv(input_data, header=TRUE, check.names=FALSE, row.names=1)
# 准备模型输入
feature_matrix <- as.matrix(expression_data)[,1:50]  # 选择前50个基因作为特征
sample_labels <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data))  # 从行名提取样本标签
# 构建Lasso回归模型
lasso_model <- glmnet(feature_matrix, sample_labels, 
                      family = "binomial", alpha=1)  # alpha=1 表示Lasso回归
                      # 图1: Lasso系数路径图
pdf(file="1_lasso_coefficient_paths.pdf", width=6, height=5.5)
plot(lasso_model)  # 展示系数随log(lambda)变化路径
dev.off()
# 交叉验证确定最优模型
set.seed(123)  # 设置随机种子保证结果可重复
cv_model <- cv.glmnet(feature_matrix, sample_labels, 
                      family="binomial", alpha=1, nfolds=10)  # 10折交叉验证
# 图2: 交叉验证误差曲线
pdf("2_cross_validation_error.pdf", width=6, height=5)
plot(cv_model)  # 展示交叉验证误差随惩罚强度变化
dev.off()
# 提取最优模型系数
optimal_lambda <- cv_model$lambda.1se  # 1个标准误规则下的lambda
model_coef <- coef(cv_model, s=optimal_lambda)  # 提取最优lambda下的系数
# 筛选重要生物标志物
selected_features <- rownames(model_coef)[which(model_coef != 0)]  # 非零系数基因
selected_features <- selected_features[-1]  # 移除截距项(第一个元素)
# 准备特征重要性数据
feature_importance <- as.matrix(model_coef)[selected_features, ]  # 提取基因系数
sorted_idx <- order(abs(feature_importance), decreasing = FALSE)  # 按绝对值升序排序
sorted_importance <- feature_importance[sorted_idx]  # 排序后的系数
# 图3: 带方向的特征重要性
pdf("3_feature_importance_directional.pdf", width=7, height=5)
barplot(sorted_importance, horiz=TRUE, las=1,
        col=ifelse(sorted_importance > 0, "tomato", "steelblue"),  # 正负系数不同颜色
        xlab="Coefficient Value", 
        main="Biomarker Impact Direction")
dev.off()
# 图4: 特征绝对重要性
pdf("4_feature_importance_absolute.pdf", width=7, height=5)
barplot(abs(sorted_importance), horiz=TRUE, las=1,
        col = "tomato",  # 统一颜色
        xlab="Absolute Coefficient Value", 
        main="Gene Importance Ranking")
dev.off()
# 保存结果
write.table(data.frame(Biomarker=selected_features, 
                       Coefficient=feature_importance),
            "selected_biomarkers.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)
# 图5: 单个基因ROC曲线
# 定义图形的颜色
bioCol=rainbow(length(selected_features), s=0.9, v=0.9)  # 生成彩虹色用于不同基因
# 提取所有基因表达数据 (前面只用了前50个基因，这里用全部181个)
expression_data <- read.csv(input_data, header=TRUE, check.names=FALSE, row.names=1)
expression_data <- expression_data[,1:181]  # 使用全部基因
# 提取样品的分组信息(对照组和实验组)
y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(expression_data))  # 提取第三个部分
y <- ifelse(y=="Control", 0, 1)  # 转换为0/1二值 (0=对照,1=实验)
# 对每个选中的基因绘制ROC曲线
aucText=c()  # 存储AUC值文本
k=0  # 计数器
# 生成PDF文件
pdf(file="5_ROC.genes.pdf", width=9, height=9)
# 循环绘制每个基因的ROC曲线
for(x in as.vector(selected_features)){
  k=k+1
  # 绘制ROC曲线
  roc1=pROC::roc(y, as.numeric(t(expression_data)[x,]))  # 计算单个基因的ROC
  if(k==1){  # 第一个基因初始化图形
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{  # 后续基因添加到同一图形
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
# 添加图例
legend("bottomright", aucText, lwd=3, bty="n", 
       cex=0.8, col=bioCol[1:k],  # 只使用已绘制基因的颜色
       inset = c(0.05, 0))  # 图例位置微调
dev.off()
# 图6: 多基因联合模型的ROC曲线
# 提取被选基因的表达数据和样本标签
selected_data <- data.frame(
  expression_data[, selected_features[1:5], drop=FALSE],  # 只取前5个基因
  Status = factor(sample_labels)  # 添加分组标签
)
# 转换分组为0/1变量
selected_data$Status <- ifelse(selected_data$Status == "DN", 1, 0)  # DN组编码为1
# 加载pROC包
library(pROC)
# 构建包含所有选定基因的模型
full_model <- glm(Status ~ ., 
                  data = selected_data,
                  family = binomial())  # 逻辑回归模型
# 计算综合预测概率
selected_data$pred_prob <- predict(full_model, type = "response")  # 预测患病概率
# 绘制ROC曲线
roc_obj <- roc(Status ~ pred_prob, data = selected_data)
# 生成PDF
pdf('6_logistical_curve.pdf', width=6.62, height=5.25)
plot(roc_obj, print.auc=TRUE)  # 绘制多基因联合ROC曲线
# AUC值表示多基因联合预测能力
# 比较单个基因AUC（图5）和联合模型AUC可评估基因组合效果
dev.off()