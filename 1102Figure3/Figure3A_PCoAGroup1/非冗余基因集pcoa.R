# ============================================================
# PCoA 分析完整脚本：基于 KEGG Pathway L2 丰度表
# 输入文件：
#   all_KEGG.PathwayL2.raw.txt  # 丰度表，行=Pathway，列=样本
#   metadata.txt                # 分组信息，行=样本，列=分组
# 输出：
#   PCoA_KEGG.pdf / PCoA_KEGG.png
#   PERMANOVA_result.txt
# ============================================================

# 1. 加载包
library(vegan)
library(ggplot2)
library(dplyr)
#all_TPM.CAZy.raw.txt
#all_eggnog.KEGG_ko.raw.txt
# 2. 读取丰度表
abund <- read.table("all_TPM.CAZy.raw.txt",
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    check.names = FALSE,
                    comment.char = "",
                    stringsAsFactors = FALSE)

# 3. 转置：样本 × Pathway
mat <- t(as.matrix(abund))

# 4. 去除全 0 的 Pathway
mat <- mat[, colSums(mat) > 0, drop = FALSE]

# 5. 读取 metadata
sub_design <- read.table("metadata.txt",
                         header = TRUE,
                         row.names = 1,
                         sep = "\t",
                         check.names = FALSE,
                         stringsAsFactors = FALSE)

# 6. 检查样本交集
samples <- intersect(rownames(sub_design), rownames(mat))
cat("丰度表样本数：", nrow(mat), "\n")
cat("metadata 样本数：", nrow(sub_design), "\n")
cat("交集样本数：", length(samples), "\n")

if (length(samples) == 0) {
  stop("丰度表和 metadata 没有共同样本，请检查样本名是否一致！")
}

# 只保留共同样本，并保证顺序一致
mat <- mat[samples, , drop = FALSE]
sub_design <- sub_design[samples, , drop = FALSE]

# 7. 计算 Bray-Curtis 距离矩阵
beta <- vegdist(mat, method = "bray")

# 8. PCoA 分析
pcoa <- cmdscale(beta, k = 2, eig = TRUE)

# 提取坐标
points <- as.data.frame(pcoa$points)
colnames(points) <- c("PCoA1", "PCoA2")

# 提取特征值并计算解释度（只考虑正特征值）
eig <- pcoa$eig
eig_pos <- eig[eig > 0]
pc1_exp <- round(100 * eig[1] / sum(eig_pos), 2)
pc2_exp <- round(100 * eig[2] / sum(eig_pos), 2)

# 9. 合并分组信息
# 请根据 metadata 中实际的分组列名修改 group_col
group_col <- "Group"   # 如果分组列叫 Treatment，请改为 "Treatment"
if (!group_col %in% colnames(sub_design)) {
  stop(paste0("metadata 中找不到分组列 '", group_col, "'，请检查列名！"))
}

points$Group <- sub_design[[group_col]]

# 设置分组因子水平（根据你的实际分组名称修改）
points$Group <- factor(points$Group,
                       levels = c("NM", "BM", "WM"))

# 10. 设置颜色（三个分组）
col <- c("#bf812d", "#543005", "#dfc27d")

# 11. PERMANOVA 分析
permanova <- adonis2(beta ~ Group,
                     data = sub_design,
                     permutations = 999)
R2 <- permanova$R2[1]
P  <- permanova$`Pr(>F)`[1]

# 输出 PERMANOVA 结果到文件
sink("PERMANOVA_result.txt")
cat("PERMANOVA 结果\n")
cat("====================\n")
print(permanova)
cat("\nR2 =", R2, "\n")
cat("P =", P, "\n")
sink()

# 12. 绘制 PCoA 图
p <- ggplot(points, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group),
               geom = "polygon",
               alpha = 0.3,
               color = NA) +
  labs(
    x = paste0("PCoA 1 (", pc1_exp, "%)"),
    y = paste0("PCoA 2 (", pc2_exp, "%)"),
    title = paste0("PERMANOVA: R² = ", round(R2, 3),
                   ", P = ", signif(P, 3))
  ) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

# 13. 显示图形
print(p)

# 14. 保存图形
ggsave("PCoA_KEGG.pdf", p, width = 6, height = 5)
ggsave("PCoA_KEGG.png", p, width = 6, height = 5, dpi = 300)

cat("\n分析完成！\n")
cat("已生成文件：\n")
cat("  PCoA_KEGG.pdf\n")
cat("  PCoA_KEGG.png\n")
cat("  PERMANOVA_result.txt\n")