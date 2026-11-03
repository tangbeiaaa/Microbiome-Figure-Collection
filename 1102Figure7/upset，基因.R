# ===============================
# 1. 加载包
# ===============================
library(tidyverse)
library(readxl)
library(dplyr)
library(ComplexUpset)
# ===============================
# 2. 读入数据
# ===============================
dat <- read.table(
  "gene_presence_absence.Rtab",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


group <- read_excel("基因组分组信息表1.xlsx")


# ===============================
# 3. 保留所有 genome（包括 Reference）
# ===============================
keep_genomes <- group$Genome
dat <- dat %>%
  select(Gene, all_of(keep_genomes))

# ===============================
# 4. 构建 group 列表（关键）
# ===============================
group_list <- split(group$Genome, group$Group)

# ===============================
# 5. 直接生成 group-level 矩阵
# ===============================
group_mat <- data.frame(Gene = dat$Gene)

for (g in names(group_list)) {
  
  genomes <- group_list[[g]]
  
  # 每个 gene 在该 group 中是否出现（任意一个 genome 有=1）
  group_mat[[g]] <- as.integer(
    rowSums(dat[, genomes, drop = FALSE]) > 0
  )
}

# ===============================
# 6. 查看结果
# ===============================
head(group_mat)
write.csv(group_mat, "group_mat.csv", row.names = FALSE)
data_mat_sub <- group_mat[, -1]  # 假设第一列是 primary_cluster
data_mat_sub <- data_mat_sub[, c("BM", "NM", "WM")]
p <- upset(
  data_mat_sub,
  
  intersect = colnames(data_mat_sub),
  
  # ⭐关键：关闭自动排序
  sort_sets = FALSE,
  
  min_size = 0,
  width_ratio = 0.1,
  height_ratio = 0.7,
  
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      text = list(size = 3)
    )
  ),
  
  themes = NULL
)
p

wm_only <- group_mat %>%
  filter(WM == 1 & BM == 0 & NM == 0 & Reference == 0)

write.csv(wm_only, "WM_only_genes.csv", row.names = FALSE)

bm_only <- group_mat %>%
  filter(BM == 1 & WM == 0 & NM == 0 & Reference == 0)

write.csv(bm_only, "BM_only_genes.csv", row.names = FALSE)
nm_only <- group_mat %>%
  filter(NM == 1 & WM == 0 & BM == 0 & Reference == 0)

write.csv(nm_only, "NM_only_genes.csv", row.names = FALSE)

wm_bm_only <- group_mat %>%
  filter(WM == 1 & BM == 1 & NM == 0 & Reference == 0)
write.csv(wm_bm_only, "WM_BM_genes.csv", row.names = FALSE)
core_genes <- group_mat %>%
  filter(WM == 1 & BM == 1 & NM == 1)
write.csv(core_genes, "core_genes.csv", row.names = FALSE)