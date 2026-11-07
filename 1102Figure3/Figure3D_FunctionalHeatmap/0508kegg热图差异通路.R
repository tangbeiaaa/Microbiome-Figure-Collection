# ====================================================
# 0. 加载包
# ====================================================
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(rstatix)
library(pheatmap)

# ====================================================
# 1. 读取数据
# ====================================================
dat <- read.delim("all_TPM.CAZyL1.raw.spf",
                  header = TRUE,
                  check.names = FALSE)

rownames(dat) <- dat$CAZy
dat <- dat[, -1, drop = FALSE]

# ====================================================
# 2. log2 转换
# ====================================================
dat_log <- log2(dat + 1)


########

# ====================================================
# 3. 读取分组
# ====================================================
group_info <- read.delim("metadata.txt",
                         header = TRUE,
                         stringsAsFactors = FALSE)

desired_order <- c("WM", "NM", "BM")
group_info$Group <- factor(group_info$Group, levels = desired_order)

# ====================================================
# 4. 转长表
# ====================================================
long_dat <- dat_log %>%
  as.data.frame() %>%
  mutate(CAZy = rownames(.)) %>%
  pivot_longer(-CAZy,
               names_to = "Sample",
               values_to = "Abundance") %>%
  left_join(group_info,
            by = c("Sample" = "Samples")) %>%
  filter(!is.na(Group))

# ====================================================
# 5. Kruskal-Wallis + FDR
# ====================================================
kw_res <- long_dat %>%
  group_by(CAZy) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"),
         sig = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ ""
         ))

# ====================================================
# 6. pairwise Wilcoxon + FDR（重点输出文件）
# ====================================================
pairwise_res <- long_dat %>%
  group_by(CAZy) %>%
  pairwise_wilcox_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  mutate(sig = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# ====================================================
# 7. 保存 pairwise 结果（你要的文件）
# ====================================================
write.csv(pairwise_res,
          "CAZy_pairwise_Wilcoxon_FDR_results.csv",
          row.names = FALSE)

# ====================================================
# 8. 构建 heatmap 数据
# ====================================================
mat_abun <- long_dat %>%
  group_by(CAZy, Group) %>%
  summarise(mean = mean(Abundance, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Group,
              values_from = mean) %>%
  column_to_rownames("CAZy") %>%
  as.matrix()

# ====================================================
# 9. Z-score
# ====================================================
mat_z <- t(scale(t(mat_abun)))
mat_z[is.na(mat_z)] <- 0

# ====================================================
# 10. heatmap 行标签（KW星号）
# ====================================================
row_labels <- rownames(mat_z)

row_labels <- paste0(
  row_labels,
  " ",
  kw_res$sig[match(rownames(mat_z), kw_res$CAZy)]
)

# ====================================================
# 11. 画 heatmap
# ====================================================
pheatmap(
  mat_z,
  color = colorRampPalette(
    c("#71B0FD", "white", "#Fd9c99")
  )(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  labels_row = row_labels,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  main = "CAZy"
)

#####kegg
# ====================================================
# 0. packages
# ====================================================
library(dplyr)
library(tidyr)
library(tibble)
library(rstatix)
library(pheatmap)

# ====================================================
# 1. read data
# ====================================================
dat <- read.delim("all_KEGG.PathwayL2.raw.txt",
                  header = TRUE,
                  check.names = FALSE)

rownames(dat) <- dat[[1]]
dat <- dat[, -1, drop = FALSE]

# log transform
dat_log <- log2(dat + 1)

# ====================================================
# 2. group info
# ====================================================
group_info <- read.delim("metadata.txt",
                         header = TRUE,
                         stringsAsFactors = FALSE)

desired_order <- c("WM", "NM", "BM")
group_info$Group <- factor(group_info$Group, levels = desired_order)

# ====================================================
# 3. long format
# ====================================================
long_dat <- dat_log %>%
  as.data.frame() %>%
  mutate(Pathway = rownames(.)) %>%
  pivot_longer(-Pathway,
               names_to = "Sample",
               values_to = "Abundance") %>%
  left_join(group_info,
            by = c("Sample" = "Samples")) %>%
  filter(!is.na(Group),
         !is.na(Abundance))

# ====================================================
# 4. ⭐ 自动过滤无效 Pathway（关键）
# ====================================================
long_dat <- long_dat %>%
  group_by(Pathway) %>%
  filter(
    n_distinct(Group[Abundance > 0]) >= 2   # 至少2组有信号
  ) %>%
  ungroup()

# ====================================================
# 5. Kruskal-Wallis + FDR
# ====================================================
kw_res <- long_dat %>%
  group_by(Pathway) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"),
         sig = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ ""
         ))

# ====================================================
# 6. pairwise Wilcoxon + FDR
# ====================================================
pairwise_res <- long_dat %>%
  group_by(Pathway) %>%
  pairwise_wilcox_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  mutate(sig = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# save table
write.csv(pairwise_res,
          "KEGG_Pathway_pairwise_Wilcoxon_FDR.csv",
          row.names = FALSE)

# ====================================================
# 7. mean matrix
# ====================================================
mat_abun <- long_dat %>%
  group_by(Pathway, Group) %>%
  summarise(mean = mean(Abundance, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Group,
              values_from = mean) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

# ====================================================
# 8. ⭐ Top15（按平均丰度）
# ====================================================
top_pathways <- order(rowMeans(mat_abun, na.rm = TRUE),
                      decreasing = TRUE)[1:min(15, nrow(mat_abun))]

mat_abun <- mat_abun[top_pathways, , drop = FALSE]

# ====================================================
# 9. Z-score
# ====================================================
mat_z <- t(scale(t(mat_abun)))
mat_z[is.na(mat_z)] <- 0

# ====================================================
# 10. KW label
# ====================================================
kw_small <- kw_res %>%
  filter(Pathway %in% rownames(mat_z))

row_labels <- paste0(
  rownames(mat_z),
  " ",
  kw_small$sig[match(rownames(mat_z), kw_small$Pathway)]
)

# ====================================================
# 11. heatmap
# ====================================================
pheatmap(
  mat_z,
  color = colorRampPalette(
    c("#71B0FD", "white", "#Fd9c99")
  )(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  labels_row = row_labels,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top15 KEGG Pathways"
)




######Lipid metabolism
library(dplyr)
library(tidyr)
library(tibble)
library(rstatix)
library(pheatmap)

# ====================================================
# 0. 工具函数：去掉 [PATH:koXXXXX] 后缀 + 清空白
# ====================================================
strip_path <- function(x) {
  x <- sub("\\s*\\[PATH:.*\\]$", "", x)
  trimws(x)
}

# ====================================================
# 1. 读索引，提取脂质代谢 level 3 通路名（原样 + 去括号两版）
# ====================================================
ko_index <- read.delim("KO1-4.txt",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       check.names = FALSE)

lipid_idx_full <- ko_index %>%
  filter(PathwayL1 == "Metabolism",
         grepl("Lipid metabolism", PathwayL2, ignore.case = TRUE)) %>%
  pull(Pathway) %>%
  unique()
lipid_idx_short <- strip_path(lipid_idx_full)

cat("索引中脂质代谢 level 3 通路数：", length(lipid_idx_full), "\n")
print(head(lipid_idx_full, 5))

if (length(lipid_idx_full) == 0) {
  stop("索引里没匹配到 Lipid metabolism，请检查 PathwayL1/PathwayL2 写法")
}

# ====================================================
# 2. 读 level 3 丰度表
# ====================================================
dat <- read.delim("all_KEGG.Pathway.raw.txt",
                  header = TRUE,
                  check.names = FALSE)
rownames(dat) <- dat[[1]]
dat <- dat[, -1, drop = FALSE]

main_full  <- rownames(dat)         # 原样行名
main_short <- strip_path(main_full) # 去括号行名

# ====================================================
# 3. 匹配：先原样，失败再去括号
# ====================================================
keep <- main_full %in% lipid_idx_full
match_mode <- "full"
if (sum(keep) == 0) {
  keep <- main_short %in% lipid_idx_short
  match_mode <- "short"
}

cat("匹配模式：", match_mode, "｜匹配行数：", sum(keep), "\n")

if (sum(keep) == 0) {
  cat("\n主数据行名示例：\n"); print(head(main_full, 8))
  cat("\n索引脂质通路名示例：\n"); print(head(lipid_idx_full, 8))
  stop("主数据没有任何行匹配脂质代谢通路，请对照上面示例核对命名差异")
}

dat_lipid <- dat[keep, , drop = FALSE]
rownames(dat_lipid) <- main_short[keep]   # 图上统一显示短名

# ====================================================
# 4. log2 转换 + 分组信息
# ====================================================
dat_log <- log2(dat_lipid + 1)

group_info <- read.delim("metadata.txt",
                         header = TRUE,
                         stringsAsFactors = FALSE)
desired_order <- c("WM", "NM", "BM")
group_info$Group <- factor(group_info$Group, levels = desired_order)

# ====================================================
# 5. 长表
# ====================================================
long_dat <- dat_log %>%
  as.data.frame() %>%
  mutate(Pathway = rownames(.)) %>%
  pivot_longer(-Pathway,
               names_to = "Sample",
               values_to = "Abundance") %>%
  left_join(group_info,
            by = c("Sample" = "Samples")) %>%
  filter(!is.na(Group), !is.na(Abundance))

# ====================================================
# 6. 自动过滤无效通路（至少 2 组有信号）
# ====================================================
long_dat <- long_dat %>%
  group_by(Pathway) %>%
  filter(n_distinct(Group[Abundance > 0]) >= 2) %>%
  ungroup()

cat("过滤后有效脂质通路数：", n_distinct(long_dat$Pathway), "\n")

# ====================================================
# 7. Kruskal-Wallis + FDR
# ====================================================
kw_res <- long_dat %>%
  group_by(Pathway) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"),
         sig = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ ""
         ))

# ====================================================
# 8. 两两 Wilcoxon + FDR
# ====================================================
pairwise_res <- long_dat %>%
  group_by(Pathway) %>%
  pairwise_wilcox_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  mutate(sig = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

write.csv(kw_res,       "Lipid_KEGG_KruskalWallis_FDR.csv", row.names = FALSE)
write.csv(pairwise_res, "Lipid_KEGG_pairwise_Wilcoxon_FDR.csv", row.names = FALSE)

# ====================================================
# 9. 组均值矩阵（列顺序 WM / NM / BM）
# ====================================================
mat_abun <- long_dat %>%
  group_by(Pathway, Group) %>%
  summarise(mean = mean(Abundance, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean)

present_groups <- intersect(desired_order, colnames(mat_abun))
mat_abun <- mat_abun %>%
  column_to_rownames("Pathway") %>%
  as.matrix()
mat_abun <- mat_abun[, present_groups, drop = FALSE]
mat_abun[is.na(mat_abun)] <- 0

# ====================================================
# 10. Top N（按平均丰度）
# ====================================================
n_show <- min(30, nrow(mat_abun))
top_pathways <- order(rowMeans(mat_abun, na.rm = TRUE),
                      decreasing = TRUE)[1:n_show]
mat_abun <- mat_abun[top_pathways, , drop = FALSE]

# ====================================================
# 11. Z-score
# ====================================================
mat_z <- t(scale(t(mat_abun)))
mat_z[is.na(mat_z)] <- 0

# ====================================================
# 12. 行标签加 KW 显著性
# ====================================================
kw_small <- kw_res %>%
  filter(Pathway %in% rownames(mat_z))

row_labels <- paste0(
  rownames(mat_z),
  " ",
  kw_small$sig[match(rownames(mat_z), kw_small$Pathway)]
)

# ====================================================
# 13. 热图
# ====================================================
pheatmap(
  mat_z,
  color = colorRampPalette(
    c("#71B0FD", "white", "#Fd9c99")
  )(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  labels_row = row_labels,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Lipid Metabolism KEGG Pathways (level 3, Top by mean abundance)"
)