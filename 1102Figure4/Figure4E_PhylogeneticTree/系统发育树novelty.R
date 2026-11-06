# =========================
# 加载包
# =========================
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(picante)
library(ggnewscale)
library(grid)
library(cowplot)
library(readxl)

# =========================
# 读取数据
# =========================
CRBD_metadata <- read_excel("MAG_tax.xlsx")
CRBD_tree <- read.tree("tax.unrootedAll.tree")
color_code <- read_excel("color_code.xlsx")

# =========================
# 数据整理
# =========================
CRBD_Species_metadata <- CRBD_metadata %>%
  mutate(
    novelty = ifelse(is.na(novelty), "known", novelty),
    label = ID,
    species_size = cluster_members,
    taxonomy = Phylum,
    taxonomy = gsub("d__|p__|c__|o__|f__|g__|s__", "", taxonomy),
    taxonomy = trimws(taxonomy)
  )

# taxonomy top12 + Others
Tax_freq <- CRBD_Species_metadata %>%
  count(taxonomy) %>%
  arrange(desc(n))

top12 <- head(Tax_freq$taxonomy, 12)

CRBD_Species_metadata <- CRBD_Species_metadata %>%
  mutate(
    taxonomy = ifelse(taxonomy %in% top12, taxonomy, "Others"),
    taxonomy = factor(taxonomy, levels = c(top12, "Others"))
  )

# =========================
# 绘图数据（不改结构）
# =========================
description_plot <- CRBD_Species_metadata %>%
  transmute(label, taxonomy_bar = taxonomy, count = 1)

species_size_plot <- CRBD_Species_metadata %>%
  transmute(
    label,
    log10size = log10(species_size)
  )

novelty_plot <- CRBD_Species_metadata %>%
  transmute(label, novelty_bar = novelty, count = 1)

# =========================
# 计算 species_size 的真实最小值和最大值（用于图例）
# =========================
min_size <- min(CRBD_Species_metadata$species_size, na.rm = TRUE)
max_size <- max(CRBD_Species_metadata$species_size, na.rm = TRUE)

# =========================
# 颜色
# =========================
tax_color <- color_code$Color
names(tax_color) <- color_code$Label
tax_color["Others"] <- "#CCCCCC"

# =========================
# tree (修改比例尺宽度为 0.4)
# =========================
CRBD_Species_metadata_for_tree <- CRBD_Species_metadata %>%
  column_to_rownames("label")

p1 <- ggtree(CRBD_tree, layout = 'fan', right = TRUE, open.angle = 15) %<+%
  CRBD_Species_metadata_for_tree

p1$data$Source[is.na(p1$data$Source)] <- "Public datasets"

p1 <- p1 +
  geom_aline(aes(color = Source), linetype = 'solid', size = 0.2, key_glyph = 'rect') +
  scale_color_manual(values = c(
    'Study samples' = '#61BFBE',
    'Public datasets' = '#F0CEA0'
  ), name = "Source") +
  geom_treescale(x = 0.6, y = -0.6, width = 0.4, color = "black", fontsize = 3) +   # 宽度改为 0.4
  theme(
    legend.position = 'left',
    legend.key = element_rect(fill = 'white', color = NA),
    legend.key.spacing.y = unit(1.5, 'mm'),
    legend.key.size = unit(1, 'line')
  )

# =========================
# taxonomy layer
# =========================
p2 <- rotate_tree(p1, angle = 170) +
  geom_fruit(
    data = description_plot,
    geom = geom_bar,
    stat = 'identity',
    width = 0.7,
    aes(y = label, x = count, fill = taxonomy_bar),
    pwidth = 0.05,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = tax_color, name = "Taxonomy")

# =========================
# size 图层 + novelty 图层
# =========================
p3 <- p2 +
  geom_fruit(
    data = species_size_plot,
    geom = geom_bar,
    stat = 'identity',
    width = 0.7,
    aes(y = label, x = log10size),
    fill = "black",
    pwidth = 0.25,
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(order = 2)) +   # 只保留 color 图例顺序，不干扰 fill
  new_scale_fill() +
  geom_fruit(
    data = novelty_plot,
    geom = geom_bar,
    stat = 'identity',
    width = 0.7,
    aes(y = label, x = count, fill = novelty_bar),
    pwidth = 0.05,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Novelty",
    values = c('novel' = '#E09220', 'known' = '#12AA9C')
  )

# =========================
# 提取右侧图例（Source, Taxonomy, Novelty）
# =========================
legend <- get_legend(
  p3 + theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )
)

# =========================
# 手动构建比例尺图例 (显示 0.4)
# =========================
tree_scale_grob <- grobTree(
  textGrob("Tree scale:", x = 0.1, y = 0.8, hjust = 0, 
           gp = gpar(fontsize = 10, fontface = "bold")),
  segmentsGrob(x0 = 0.1, x1 = 0.5, y0 = 0.5, y1 = 0.5, 
               gp = gpar(lwd = 1.5)),
  segmentsGrob(x0 = 0.1, x1 = 0.1, y0 = 0.45, y1 = 0.55),
  segmentsGrob(x0 = 0.5, x1 = 0.5, y0 = 0.45, y1 = 0.55),
  textGrob("0.4", x = 0.3, y = 0.3, gp = gpar(fontsize = 9))
)

# =========================
# 手动构建 species size 图例 (动态最小最大值，但保持 log10 形式)
# =========================
size_legend_grob <- grobTree(
  textGrob(expression(Log[10] ~ "(species size)"), 
           x = 0.1, y = 0.8, hjust = 0, 
           gp = gpar(fontsize = 10, fontface = "bold")),
  rectGrob(x = 0.15, y = 0.55, width = 0.08, height = 0.08, 
           gp = gpar(fill = "black", col = NA)),
  textGrob(paste0("Min. log10 (", min_size, "), max. log10 (", max_size, ")"), 
           x = 0.28, y = 0.55, hjust = 0, gp = gpar(fontsize = 9))
)

# =========================
# 整合最终图形
# =========================
final_plot <- ggdraw() +
  draw_plot(
    p3 + theme(legend.position = "none"),
    x = 0, y = 0, width = 0.85, height = 1
  ) +
  draw_plot(
    legend,
    x = 0.85, y = 0, width = 0.15, height = 1
  ) +
  draw_plot(
    tree_scale_grob,
    x = 0.01, y = 0.01, width = 0.2, height = 0.12
  ) +
  draw_plot(
    size_legend_grob,
    x = 0.01, y = 0.14, width = 0.25, height = 0.12
  )

# 显示图形
print(final_plot)
# 保存（可选）
# ggsave("tree_with_manual_legends.pdf", final_plot, width = 12, height = 9, dpi = 300)