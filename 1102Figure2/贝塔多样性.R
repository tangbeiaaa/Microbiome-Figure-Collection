# =========================
# PCoA 分析及可视化
# =========================

library(vegan)      # α/β 多样性分析
library(readxl)     # 读取 Excel
library(readr)      # 读取 CSV/TXT
library(reshape2)   # 转矩阵
library(ggplot2)    # 可视化
library(dplyr)      # 数据处理
library(ggpubr)     # 添加统计显著性
# -------------------------
# 1. 读取距离矩阵
# -------------------------
# 1. 读取表格，第一列是 MAG 名或其他标识

读取距离矩阵

beta <- read.table("bray_curtis.txt", 
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t", 
                   comment.char = "")

#读取样本表型数据
sub_design<-read.table("metadata.txt",header=T,row.names=1,sep="\t")
# 转成 data.frame
sub_design <- as.data.frame(sub_design)
# 查看结构
str(sub_design)

head(sub_design)
str(sub_design)
dim(sub_design)
colnames(sub_design)

# -------------------------
# 2. 取样本设计表与距离矩阵交集
# -------------------------
# 取 sub_design 中存在于 beta 的样本
samples <- rownames(sub_design)
samples <- samples[samples %in% rownames(beta)]

# 构建子集 beta 矩阵
sub_beta <- beta[samples, samples]

# 查看维度确认
dim(sub_beta)
head(sub_beta[,1:6])
#idx <- rownames(sub_design) %in% rownames(beta)
#sub_design <- sub_design[idx, ]
#sub_beta <- beta[rownames(sub_design), rownames(sub_design)]

# -------------------------
# 3. PCoA分析（基于距离矩阵）
# -------------------------
# k = 3：输出前三个主坐标
# eig = TRUE：返回特征值
pcoa <- cmdscale(sub_beta, k = 3, eig = TRUE)
str(pcoa)

# 提取前三个主坐标
points <- as.data.frame(pcoa$points)
colnames(points) <- c("x", "y", "z")
head(points)

# 提取特征值
eig <- pcoa$eig

# -------------------------
# 4. 合并样本表型信息
# -------------------------
points$Group <- sub_design$Group[match(rownames(points), rownames(sub_design))]
head(points)

# 设置因子
points$Group <- factor(points$Group, 
                       levels = c("NM", "BM", "WM"))
points$Treatment <- factor(points$Group, 
                           levels = c("NM", "BM", "WM"))

# -------------------------
# 5. 设置绘图颜色
# -------------------------
col <- c("#bf812d","#543005","#dfc27d")
# -------------------------
# 6. 绘制PCoA散点图
# -------------------------
p <- ggplot(points, aes(x = x, y = y, color = Group, shape =Group)) + 
  geom_point(size = 3) +
  labs(
    x = paste0("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)"),
    y = paste0("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)")
  ) +
  scale_colour_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )

# -------------------------
# 7. 展示图形
# -------------------------
print(p)


p <- ggplot(points, aes(x = x, y = y, color = Group, shape = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +  # 半透明底色椭圆
  labs(
    x = paste0("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)"),
    y = paste0("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)")
  ) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

print(p)


#
# 7. PERMANOVA分析
# -------------------------
dune_div <- adonis2(sub_beta ~ Group, data = sub_design, permutations = 999, method = "bray")

# 输出 PERMANOVA 结果
print(dune_div)

# -------------------------
# 8. 绘制带统计显著性的PCoA图
# -------------------------
R2 <- dune_div$R2[1]
P  <- dune_div$`Pr(>F)`[1]

p <- ggplot(points, aes(x = x, y = y, color = Group, shape = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group),
               geom = "polygon",
               alpha = 0.3,
               color = NA) +
  labs(
    x = paste0("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)"),
    y = paste0("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)"),
    title = paste0(
      "PERMANOVA: R² = ",
      round(R2, 3),
      ", P = ",
      signif(P, 3)
    )
  ) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  theme_classic()

print(p)