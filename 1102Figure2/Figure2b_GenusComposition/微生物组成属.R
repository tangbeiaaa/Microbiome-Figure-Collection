library(tidyverse)
library(patchwork)
library(scales)
library(RColorBrewer)

# ================= 1. Read abundance table =================
species <- read.delim(
  "bracken.G.0.2",
  row.names = 1,
  check.names = FALSE
)

# 【关键修复1】把行名里的 "g__" 前缀去掉，防止跟颜色对不上
rownames(species) <- sub("^g__", "", rownames(species))

# ================= 2. Read metadata =================
metadata <- read.delim(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

group_df <- metadata %>%
  select(Samples, Group) %>%
  rename(Sample = Samples)

group_df <- group_df %>%
  filter(Sample %in% colnames(species))

group_df$Group <- factor(group_df$Group, levels = c("WM", "NM", "BM"))

# ================= 3. Select Top10 taxa =================
taxa_clean <- species %>%
  rownames_to_column("Taxa") %>%
  mutate(total = rowSums(across(-Taxa))) %>%
  arrange(desc(total)) %>%
  slice(1:10) %>%
  select(-total)

# ================= 4. Merge low abundance taxa =================
others <- species %>%
  filter(!rownames(.) %in% taxa_clean$Taxa) %>%
  summarise(across(everything(), sum)) %>%
  mutate(Taxa = "Others")

# ================= 5. Calculate relative abundance =================
df_taxa <- bind_rows(taxa_clean, others) %>%
  mutate(across(-Taxa, ~ .x / sum(.x))) %>%
  pivot_longer(
    -Taxa,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(group_df, by = "Sample")

# ================= 6. Sample order =================
# 【关键修复2】动态获取丰度第一的属作为排序依据，而不是写死为门
target_taxa <- taxa_clean$Taxa[1]
cat("用于排序的 Top1 属是：", target_taxa, "\n")

sample_order_df <- df_taxa %>%
  filter(Taxa == target_taxa) %>%
  arrange(Group, desc(Abundance)) %>%
  pull(Sample)

df_taxa$Sample <- factor(
  df_taxa$Sample,
  levels = sample_order_df
)

# ================= 7. Taxa colors =================
color_map <- c(
  "Corynebacterium" = "#3288bd",
  "Trueperella"     = "#66c2a5",
  "Pseudomonas"     = "#abdda4",
  "Proteus"         = "#fdae61",
  "Mammaliicoccus"  = "#9e0142",
  "Enterococcus"    = "#d53e4f",
  "Streptomyces"    = "#E9792F",
  "Aerococcus"      = "#F6C8A8",
  "Arthrobacter"    = "#F9CA34",
  "Bacteroides"     = "#f46d43",
  "Others"          = "#ffffbf"
)

# 【关键修复3】确保所有实际存在的属都有颜色，如果没有自动补全
all_taxa <- c(taxa_clean$Taxa, "Others")
missing_taxa <- setdiff(all_taxa, names(color_map))
if (length(missing_taxa) > 0) {
  cat("警告：以下属没有指定颜色，已自动补全：\n")
  print(missing_taxa)
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(missing_taxa))
  names(extra_colors) <- missing_taxa
  color_map <- c(color_map, extra_colors)
}
# 保证顺序与数据完全一致
# 按全样本平均相对丰度降序排列，但 Others 强制放最后
mean_order <- df_taxa %>%
  group_by(Taxa) %>%
  summarise(Mean = mean(Abundance), .groups = "drop") %>%
  arrange(desc(Mean)) %>%
  pull(Taxa)

# 把 Others 挪到末尾（若存在）
mean_order <- c(setdiff(mean_order, "Others"), intersect(mean_order, "Others"))

# 安全检查
if (!setequal(mean_order, names(color_map))) {
  stop("color_map 与数据 Taxa 不一致，缺少：",
       paste(setdiff(mean_order, names(color_map)), collapse = ", "),
       "；多余：",
       paste(setdiff(names(color_map), mean_order), collapse = ", "))
}

color_map <- color_map[mean_order]
# ================= 8. Group colors =================
Treatments_colors <- c(
  "WM" = "#dfc27d",
  "NM" = "#bf812d",
  "BM" = "#543005"
)

# ================= 9. Top annotation bar =================
plot_top_bar <- function(group_name, color_map_treat) {
  
  n_samples <- sum(group_df$Group == group_name)
  
  bar_data <- data.frame(
    x = 1:n_samples,
    y = 1,
    Group = factor(group_name, levels = names(color_map_treat))
  )
  
  ggplot(bar_data, aes(x = x, y = y, fill = Group)) +
    geom_tile() +
    scale_fill_manual(values = color_map_treat) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.margin = margin(t = 5, r = 0, b = 5, l = 0),
      legend.position = "none"
    )
}

# ================= 10. Area plot =================
plot_musk_area <- function(group_name, data, color_vec, position = "middle") {
  
  sub_data <- data %>% filter(Group == group_name)
  
  sub_data$Taxa <- factor(sub_data$Taxa, levels = names(color_vec))
  sub_data$x_numeric <- as.numeric(sub_data$Sample)
  
  p <- ggplot(sub_data, aes(x = x_numeric, y = Abundance, fill = Taxa, group = Taxa)) +
    geom_area(
      position = "fill",
      alpha = 0.85,
      colour = NA
    ) +
    scale_fill_manual(values = color_vec, name = "Genus") +
    scale_y_continuous(expand = c(0, 0), labels = percent) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    ) +
    labs(x = NULL, y = NULL)
  
  if(position == "left") {
    p <- p +
      labs(y = "Relative Abundance") +
      theme(axis.line.y = element_line(colour = "black"))
  } else {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  }
  
  return(p)
}

# ================= 11. Generate plots =================
t_WM <- plot_top_bar("WM", Treatments_colors)
a_WM <- plot_musk_area("WM", df_taxa, color_map, position = "left")

t_NM <- plot_top_bar("NM", Treatments_colors)
a_NM <- plot_musk_area("NM", df_taxa, color_map, position = "middle")

t_BM <- plot_top_bar("BM", Treatments_colors)
a_BM <- plot_musk_area("BM", df_taxa, color_map, position = "middle")

# ================= 12. Combine plots =================
col_WM <- t_WM / a_WM + plot_layout(heights = c(0.4, 20))
col_NM <- t_NM / a_NM + plot_layout(heights = c(0.4, 20))
col_BM <- t_BM / a_BM + plot_layout(heights = c(0.4, 20))

final_plot <- (col_WM | col_NM | col_BM) +
  plot_layout(guides = "collect", widths = c(10, 9, 11)) &
  theme(
    legend.position = "right",
    panel.spacing = unit(0, "cm")
  )

# ================= 13. Show plot =================
print(final_plot)

# ================= 14. Export PDF =================
ggsave(
  "Stacked_Bar_Genus_NoGap.pdf",
  final_plot,
  width = 12,
  height = 6,
  device = grDevices::cairo_pdf
)

# ================= 15. Export group mean relative abundance =================
group_mean <- df_taxa %>%
  group_by(Group, Taxa) %>%
  summarise(
    Mean = mean(Abundance),
    SD = sd(Abundance),
    .groups = "drop"
  ) %>%
  mutate(
    Mean = Mean * 100,
    SD = SD * 100
  )

write.csv(group_mean, "Group_mean_relative_abundance_long.csv", row.names = FALSE, quote = FALSE)

group_mean_wide <- group_mean %>%
  select(Taxa, Group, Mean) %>%
  pivot_wider(names_from = Group, values_from = Mean)

write.csv(group_mean_wide, "Group_mean_relative_abundance_wide.csv", row.names = FALSE, quote = FALSE)

print(group_mean_wide)

# ================= 16. Overall relative abundance =================
overall_mean <- df_taxa %>%
  group_by(Taxa) %>%
  summarise(
    Mean = mean(Abundance),
    SD = sd(Abundance),
    .groups = "drop"
  ) %>%
  mutate(
    Mean = Mean * 100,
    SD = SD * 100
  ) %>%
  arrange(desc(Mean))

write.csv(overall_mean, "Overall_relative_abundance.csv", row.names = FALSE, quote = FALSE)

print(overall_mean)