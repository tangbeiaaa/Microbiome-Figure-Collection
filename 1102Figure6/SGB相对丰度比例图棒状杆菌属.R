
library(ggplot2)
library(dplyr)
library(cowplot)
library(readxl)
library(tidyr)

# ------------------ 1. 读取数据 ------------------
abund <- read.delim(
  "group_mean_mag.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

anno <- read_excel("rename_cor.xlsx")

# 检查必要列
if (!all(c("OTUID", "bad_musk", "good_musk", "white_musk") %in% names(abund))) {
  stop("group_mean_mag.txt 缺少必要的丰度列。")
}

if (!all(c("ID", "SGB") %in% names(anno))) {
  stop("rename_cor.xlsx 必须包含 ID 和 SGB 两列。")
}

# ------------------ 2. 统一 MAG ID 格式 ------------------
normalize_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.fa$", "", x)
}

abund$ID_clean <- normalize_id(abund$OTUID)
anno$ID_clean <- normalize_id(anno$ID)

# 去掉空值和 unmapped
anno_clean <- anno %>%
  filter(
    !is.na(SGB),
    SGB != "",
    !is.na(ID_clean),
    ID_clean != ""
  ) %>%
  distinct(ID_clean, .keep_all = TRUE)

abund_clean <- abund %>%
  filter(
    !is.na(ID_clean),
    ID_clean != "",
    ID_clean != "unmapped"
  ) %>%
  mutate(
    across(
      c(bad_musk, good_musk, white_musk),
      as.numeric
    )
  )

# ------------------ 3. 只保留 rename_cor 中的 MAG ------------------
merged_data <- abund_clean %>%
  inner_join(
    anno_clean %>% select(ID_clean, SGB),
    by = "ID_clean"
  )

cat("丰度表中的 MAG 数量：", nrow(abund_clean), "\n")
cat("成功匹配注释的 MAG 数量：", nrow(merged_data), "\n")
cat("匹配到的 SGB 数量：", n_distinct(merged_data$SGB), "\n")

if (nrow(merged_data) == 0) {
  stop("没有匹配到 MAG，请检查两个文件中的 ID 是否一致。")
}

# ------------------ 4. 相同 SGB 的丰度求和 ------------------

# ------------------ 4. 相同 SGB 的 MAG 丰度求平均 ------------------
sgb_abund <- merged_data %>%
  group_by(SGB) %>%
  summarise(
    bad_musk = mean(bad_musk, na.rm = TRUE),
    good_musk = mean(good_musk, na.rm = TRUE),
    white_musk = mean(white_musk, na.rm = TRUE),
    MAG_count = n(),
    .groups = "drop"
  )


# ------------------ 5. 按三组平均丰度排序 ------------------
sgb_abund <- sgb_abund %>%
  mutate(
    Mean_Abundance = rowMeans(
      across(c(bad_musk, good_musk, white_musk)),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(Mean_Abundance))

# 取丰度最高的 20 个 SGB
top20_sgb <- sgb_abund %>%
  slice_head(n = 20)

print(top20_sgb)

# 最高丰度的 SGB 显示在最上方
sgb_order <- rev(top20_sgb$SGB)

# ------------------ 6. 整理成长格式 ------------------
plot_data <- top20_sgb %>%
  select(SGB, bad_musk, good_musk, white_musk) %>%
  pivot_longer(
    cols = c(bad_musk, good_musk, white_musk),
    names_to = "Group",
    values_to = "Abundance"
  )

plot_data$SGB <- factor(
  plot_data$SGB,
  levels = sgb_order
)

# ------------------ 7. 分组名称与颜色 ------------------
# 假定：
# bad_musk   = BM
# good_musk  = NM
# white_musk = WM

group_map <- c(
  bad_musk = "BM",
  good_musk = "NM",
  white_musk = "WM"
)

plot_data$Group <- factor(
  plot_data$Group,
  levels = c("white_musk", "good_musk", "bad_musk"),
  labels = c("WM", "NM", "BM")
)

colors <- c(
  "WM" = "#F0E4B8",
  "NM" = "#C8963E",
  "BM" = "#4A2E14"
)

# ------------------ 8. 绘制三个子图 ------------------
groups <- c("WM", "NM", "BM")

# 🔑 关键：所有面板共用一个 y 轴上限，保证可比
y_max <- max(plot_data$Abundance, na.rm = TRUE) * 1.05

plots <- lapply(seq_along(groups), function(i) {
  
  g <- groups[i]
  dat <- plot_data %>%
    filter(Group == g)
  
  p <- ggplot(
    dat,
    aes(x = SGB, y = Abundance)
  ) +
    geom_col(
      fill = colors[g],
      width = 0.75
    ) +
    coord_flip() +
    labs(
      title = g,
      x = NULL,
      y = "relation abundance"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "sans", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10),
      plot.title = element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      ),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    scale_y_continuous(
      expand = c(0, 0),            # 去掉默认扩展，避免和 limits 冲突
      limits = c(0, y_max)         # 🔑 三组共用同一范围
    )
  
  # 仅在第一张图显示 SGB 名称
  if (i == 1) {
    p <- p +
      theme(
        axis.text.y = element_text(
          size = 9,
          color = "black",
          hjust = 1
        )
      )
  }
  
  p
})

# ------------------ 9. 拼图 ------------------
final_plot <- plot_grid(
  plotlist = plots,
  ncol = 3,
  rel_widths = c(1.5, 1, 1),
  align = "h",
  axis = "tb"
)

print(final_plot)
###################
# ------------------ 8. 合并为一张分组柱状图 ------------------
p_combined <- ggplot(
  plot_data,
  aes(x = SGB, y = Abundance, fill = Group)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = NA
  ) +
  coord_flip() +
  scale_fill_manual(
    values = colors,
    name   = "Group",
    breaks = c("WM", "NM", "BM")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "relative abundance"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "sans", size = 12),
    axis.text.y      = element_text(size = 9, color = "black", hjust = 1),
    axis.text.x      = element_text(size = 9, color = "black"),
    axis.title.x     = element_text(size = 11),
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 10)
  )

print(p_combined)

# 保存
ggsave("SGB_top20_grouped_bar.pdf", p_combined, width = 8, height = 7)
ggsave("SGB_top20_grouped_bar.png", p_combined, width = 8, height = 7, dpi = 300)