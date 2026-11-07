library(tidyverse)
library(readxl)
library(pheatmap)

# ====================================================
# 0. 工具函数
# ====================================================
strip_path <- function(x) trimws(sub("\\s*\\[PATH:.*\\]$", "", x))

# ====================================================
# 1. 读 CAZy 丰度矩阵
# ====================================================
cazy_raw <- read.delim("sum.Level1.raw.cor.txt",
                       row.names = 1, check.names = FALSE)
cat("CAZy 原始维度：", dim(cazy_raw), "\n")
cat("值范围：", range(as.matrix(cazy_raw), na.rm = TRUE), "\n")

# ====================================================
# 2. 读重命名表与分组表
# ====================================================
rename_df <- read_excel("rename_cor.xlsx") %>% as.data.frame()
colnames(rename_df) <- c("ID", "SGB")
rename_df$ID_clean <- sub("\\.fa$", "", rename_df$ID)
rename_df$ID_clean <- sub("^(Sg_|Gp_)", "", rename_df$ID_clean)

group_df <- read_excel("cor_group.xlsx") %>% as.data.frame()
colnames(group_df) <- c("ID", "group")
group_df$ID_clean <- sub("\\.fa$", "", group_df$ID)
group_df$ID_clean <- sub("^(Sg_|Gp_)", "", group_df$ID_clean)

colnames(cazy_raw) <- sub("^(Sg_|Gp_)", "", colnames(cazy_raw))

sgb_vec <- rename_df$SGB[match(colnames(cazy_raw), rename_df$ID_clean)]
sgb_vec[is.na(sgb_vec)] <- "NA"
group_vec <- group_df$group[match(colnames(cazy_raw), group_df$ID_clean)]
group_vec[is.na(group_vec)] <- "NA"

# ====================================================
# 3. 只保留指定 4 个 SGB；行不做任何过滤
# ====================================================
target_sgbs <- c("SGB90_1", "SGB85_0", "SGB84_2", "SGB83_1")
keep_idx <- sgb_vec %in% target_sgbs

cazy_mat  <- cazy_raw[, keep_idx, drop = FALSE]
sgb_vec   <- sgb_vec[keep_idx]
group_vec <- group_vec[keep_idx]

cat("按 SGB 过滤后维度：", dim(cazy_mat), "\n")
print(table(sgb_vec))
if (ncol(cazy_mat) == 0) stop("4 个目标 SGB 一列都没匹配上")

# —— 转数值矩阵 + NA 置 0（不做任何行过滤）——
mat_cnt <- as.matrix(cazy_mat)
mat_cnt <- apply(mat_cnt, 2, function(x)
  suppressWarnings(as.numeric(as.character(x))))
rownames(mat_cnt) <- rownames(cazy_mat)
mat_cnt[is.na(mat_cnt)] <- 0

cat("保留全部 CAZy 行，行数：", nrow(mat_cnt), "\n")
rownames(mat_cnt) <- make.unique(rownames(mat_cnt))

# ====================================================
# 4. CPM 归一化 + log2 变换
# ====================================================
lib_size <- colSums(mat_cnt)
lib_size[lib_size == 0] <- 1
cpm     <- sweep(mat_cnt, 2, lib_size, "/") * 1e6
mat_log <- log2(cpm + 1)
cat("CPM+log2 后维度：", dim(mat_log), "\n")

# ====================================================
# 5. 手动 row z-score；NaN（全0/常数行）置 0
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
# 8. 绘制 CAZy 热图（scale="none"，直接输出到 Plots，不保存）
# ====================================================
n_row <- nrow(z)
n_col <- ncol(z)
cat("绘制矩阵规模：", n_row, "行 ×", n_col, "列\n")

cw     <- max(8, min(20, 600 / n_col))
ch     <- if (n_row > 150) 5 else
  if (n_row > 80)  10 else 18
fs_row <- if (n_row > 150) 6 else
  if (n_row > 80)  8 else 10
fs_col <- if (n_col > 40) 6 else
  if (n_col > 20) 8 else 10

pheatmap(
  z,
  scale             = "none",          # 🔑 已手动 z，不能再 scale="row"
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = fs_row,
  fontsize_col      = fs_col,
  border_color      = "white",
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  treeheight_row    = 20,
  cellwidth         = cw,
  cellheight        = ch,
  color             = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  breaks            = seq(-3, 3, length.out = 101),
  main = "Target SGBs CAZy (Level1, all rows, no filtering)"
  # 🔑 不传 filename —— 直接显示在 Plots 面板
)