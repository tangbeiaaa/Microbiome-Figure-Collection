# ============================================================
# Corynebacterium 环形树
#   - 内环：SGB 条带（fill 映射，图例自动方块）
#   - 中环：Species 条带（fill 映射）
#   - 外环：Group 五角星
#   - 分支：按 SGB 着色（图例隐藏）
# ============================================================

library(readxl)
library(tidyverse)
library(treeio)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
library(ggnewscale)
library(ggstar)

# ---------- 1. 读树 + 去前缀 ----------
tree <- read.tree("Corynebacterium_tree.treefile")
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- sub("^(Sg_|Gp_)", "", tree$tip.label)

# ---------- 2. 读数据 ----------
sgb_info <- read_excel("Corynebacterium_SGB.xlsx") %>% as.data.frame()
colnames(sgb_info)[1:5] <- c("genome", "SGB", "Species", "group", "res")
sgb_info$tip_label <- sub("\\.fa$", "", sgb_info$genome)
sgb_info$tip_label <- sub("^(Sg_|Gp_)", "", sgb_info$tip_label)

tip_meta <- sgb_info %>%
  select(tip_label, SGB, Species, group) %>%
  filter(tip_label %in% tree$tip.label) %>%
  distinct()

tip_meta$group[is.na(tip_meta$group) | tip_meta$group == ""] <- "NA"
tip_meta$Species[is.na(tip_meta$Species) | tip_meta$Species == ""] <- "NA"
tip_meta$SGB[is.na(tip_meta$SGB) | tip_meta$SGB == ""] <- "NA"

# ---------- 3. 配色 ----------

# SGB 配色
sgb_palette <- c(
  "#3188ba", "#65bfa4", "#3D9F3C", "#abd7a2",
  "#f5f8be", "#fbae61", "#E9E483", "#A1dce6",
  "#F4B6C2", "#BDACD6", "#A1CDEC", "#f1c8a8"
)

sgb_levels <- sort(unique(tip_meta$SGB))
n_sgb <- length(sgb_levels)
if (n_sgb <= length(sgb_palette)) {
  sgb_colors <- sgb_palette[1:n_sgb]
} else {
  sgb_colors <- colorRampPalette(sgb_palette)(n_sgb)
}
names(sgb_colors) <- sgb_levels

# Species 配色
species_levels <- unique(tip_meta$Species)
if (length(species_levels) <= 8) {
  species_colors <- brewer.pal(max(3, length(species_levels)), "Set2")
} else {
  species_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(species_levels))
}
names(species_colors) <- species_levels

# Group 配色
group_colors <- c(
  "WM" = "#F0E4B8",
  "NM" = "#C8963E",
  "BM" = "#4A2E14",
  "NA" = "grey70"
)

# ---------- 4. 基础树 + 分支按 SGB 着色 ----------
p <- ggtree(tree, layout = "circular", size = 0.8,
            branch.length = "none") %<+% tip_meta

sgb_list <- split(tip_meta$tip_label, tip_meta$SGB)
p <- groupOTU(p, sgb_list, group_name = "SGB_group")

p <- p +
  aes(color = SGB_group) +
  scale_color_manual(values = sgb_colors, na.value = "grey40",
                     name = "SGB") +
  # ★ 隐藏分支颜色图例（下面用 SGB 环的图例替代）
  guides(color = "none") +
  theme(legend.position = "right")

# ---------- 5. 计算环的半径 ----------
tree_data <- p$data
max_x <- max(tree_data$x[tree_data$isTip], na.rm = TRUE)
tip_data <- tree_data %>% filter(isTip)

sgb_offset     <- max_x * 1.03
species_offset <- max_x * 1.12
group_offset   <- max_x * 1.25
bar_width      <- max_x * 0.035

# ---------- 6. 内环：SGB 条带（图例自动方块） ----------
p <- p +
  ggnewscale::new_scale_fill() +
  geom_tile(data = tip_data,
            aes(x = sgb_offset, y = y, fill = SGB),
            width = bar_width,
            height = 0.85,
            color = NA,
            inherit.aes = FALSE) +
  scale_fill_manual(values = sgb_colors, name = "SGB")

# ---------- 7. 中环：Species 条带 ----------
p <- p +
  ggnewscale::new_scale_fill() +
  geom_tile(data = tip_data,
            aes(x = species_offset, y = y, fill = Species),
            width = bar_width,
            height = 0.85,
            color = NA,
            inherit.aes = FALSE) +
  scale_fill_manual(values = species_colors, name = "Species")

# ---------- 8. 外环：Group 五角星 ----------
p <- p +
  ggnewscale::new_scale_fill() +
  geom_star(data = tip_data,
            aes(x = group_offset, y = y, fill = group),
            size = 2.5,
            starstroke = 0.2,
            starshape = 8,
            inherit.aes = FALSE) +
  scale_fill_manual(values = group_colors, name = "Group")

# ---------- 9. 输出 ----------
print(p)

ggsave("Corynebacterium_tree_final.pdf", p, width = 14, height = 14)
ggsave("Corynebacterium_tree_final.png", p, width = 14, height = 14, dpi = 300)

cat("完成\n")