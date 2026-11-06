library(tidyverse)
library(readxl)
library(pheatmap)

# ====================================================
# 0. 工具函数：去掉 [PATH:koXXXXX] 后缀 + 清空白
# ====================================================
strip_path <- function(x) {
  x <- sub("\\s*\\[PATH:.*\\]$", "", x)
  trimws(x)
}

# ====================================================
# 1. 读 KO 索引，提取脂质代谢 level 3 通路
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

# ====================================================
# 2. 读 level 3 丰度表
# ====================================================
dat <- read.delim("KEGG.Pathway.raw.cor.txt",
                  header = TRUE, check.names = FALSE)
rownames(dat) <- dat[[1]]
dat <- dat[, -1, drop = FALSE]

main_full  <- rownames(dat)
main_short <- strip_path(main_full)

# ====================================================
# 3. 匹配脂质代谢通路行
# ====================================================
keep <- main_full %in% lipid_idx_full
match_mode <- "full"
if (sum(keep) == 0) {
  keep <- main_short %in% lipid_idx_short
  match_mode <- "short"
}
cat("匹配模式：", match_mode, "｜匹配行数：", sum(keep), "\n")
if (sum(keep) == 0) stop("主数据未匹配到任何脂质代谢通路，请检查 KO1-4.txt 列名或行名格式")

dat_lipid <- dat[keep, , drop = FALSE]
rownames(dat_lipid) <- main_short[keep]

# ====================================================
# 4. 读重命名表与分组表
# ====================================================
rename_df <- read_excel("rename_cor.xlsx") %>% as.data.frame()
colnames(rename_df) <- c("ID", "SGB")
rename_df$ID_clean <- sub("\\.fa$", "", rename_df$ID)
rename_df$ID_clean <- sub("^(Sg_|Gp_)", "", rename_df$ID_clean)

group_df <- read_excel("cor_group.xlsx") %>% as.data.frame()
colnames(group_df) <- c("ID", "group")
group_df$ID_clean <- sub("\\.fa$", "", group_df$ID)
group_df$ID_clean <- sub("^(Sg_|Gp_)", "", group_df$ID_clean)

colnames(dat_lipid) <- sub("^(Sg_|Gp_)", "", colnames(dat_lipid))

sgb_vec <- rename_df$SGB[match(colnames(dat_lipid), rename_df$ID_clean)]
sgb_vec[is.na(sgb_vec)] <- "NA"
group_vec <- group_df$group[match(colnames(dat_lipid), group_df$ID_clean)]
group_vec[is.na(group_vec)] <- "NA"

# ====================================================
# 5. 只保留指定 4 个 SGB；行不做任何过滤，保留全部 16 条
# ====================================================
target_sgbs <- c("SGB90_1", "SGB85_0", "SGB84_2", "SGB83_1")
keep_idx <- sgb_vec %in% target_sgbs

dat_lipid <- dat_lipid[, keep_idx, drop = FALSE]
sgb_vec   <- sgb_vec[keep_idx]
group_vec <- group_vec[keep_idx]

cat("按 SGB 过滤后维度：", dim(dat_lipid), "\n")
if (ncol(dat_lipid) == 0) stop("4 个目标 SGB 一列都没匹配上")

# —— 转数值矩阵 + NA 置 0（不做任何行过滤，16 条全保留）——
mat_cnt <- as.matrix(dat_lipid)
mat_cnt <- apply(mat_cnt, 2, function(x)
  suppressWarnings(as.numeric(as.character(x))))
rownames(mat_cnt) <- rownames(dat_lipid)
mat_cnt[is.na(mat_cnt)] <- 0

cat("保留全部通路行，行数：", nrow(mat_cnt), "\n")
rownames(mat_cnt) <- make.unique(rownames(mat_cnt))

# ====================================================
# 5.5 CPM 归一化 + log2 变换
# ====================================================
lib_size <- colSums(mat_cnt)
lib_size[lib_size == 0] <- 1              # 防某列全 0 除 0
cpm     <- sweep(mat_cnt, 2, lib_size, "/") * 1e6
mat_log <- log2(cpm + 1)
cat("CPM+log2 后维度：", dim(mat_log), "\n")

# ====================================================
# 5.6 手动 row z-score；NaN（全0/常数行）置 0
# ====================================================
z <- t(scale(t(mat_log)))
n_nan <- sum(is.nan(z))
if (n_nan > 0) {
  cat("零方差行数（将显示为白色，z=0）：", n_nan, "\n")
  print(rownames(mat_log)[apply(mat_log, 1, sd) == 0])
}
z[is.nan(z)] <- 0
z[z >  3] <-  3
z[z < -3] <- -3

# ====================================================
# 6. 按品质组（WM→NM→BM→NA）与 SGB 排序
# ====================================================
group_order <- c("WM", "NM", "BM", "NA")
col_order <- order(factor(group_vec, levels = group_order), sgb_vec)

z         <- z[, col_order, drop = FALSE]
sgb_vec   <- sgb_vec[col_order]
group_vec <- group_vec[col_order]

annotation_col <- data.frame(
  Group = factor(group_vec, levels = group_order),
  SGB   = factor(sgb_vec, levels = sort(unique(sgb_vec)))
)
rownames(annotation_col) <- colnames(z)

# ====================================================
# 7. 配色
# ====================================================
group_colors <- c(
  "WM" = "#F0E4B8", "NM" = "#C8963E",
  "BM" = "#4A2E14", "NA" = "grey85"
)

sgb_palette <- c("#3D9F3C", "#f5f8be", "#fbae61", "#A1dce6")
sgb_unique  <- sort(unique(sgb_vec))
sgb_colors  <- sgb_palette[seq_along(sgb_unique)]
names(sgb_colors) <- sgb_unique
if ("NA" %in% sgb_unique) sgb_colors["NA"] <- "grey85"

ann_colors <- list(Group = group_colors, SGB = sgb_colors)

# ====================================================
# 8. 绘制脂质代谢通路热图（scale = "none"，因为已手动 z）
# ====================================================
n_col  <- ncol(z)
cw     <- max(8, min(20, 600 / n_col))
fs_col <- ifelse(n_col > 40, 6, ifelse(n_col > 20, 8, 10))

pheatmap(
  z,
  scale             = "none",          # 🔑 已经手动 z 过，不能再 scale="row"
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 10,
  fontsize_col      = fs_col,
  border_color      = "white",
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  treeheight_row    = 20,
  cellwidth         = cw,
  cellheight        = 18,
  color             = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  breaks            = seq(-3, 3, length.out = 101),   # 与截断范围一致
  main = "Target SGBs Lipid Metabolism Pathways (all 16, no filtering)"
)