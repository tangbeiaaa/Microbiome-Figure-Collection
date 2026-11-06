# 清除环境
rm(list = ls())

# 加载所需包（若无则安装）
packages <- c("ggplot2", "dplyr", "RColorBrewer", "gridExtra", "cowplot")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
    library(pkg, character.only = TRUE)
  }
}

# ==================== 1. 读取数据 ====================
input_file <- "quality_report.tsv"  # 请根据实际路径修改
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df[, c("Completeness", "Contamination", "Genome_Size")]
df$Completeness <- as.numeric(df$Completeness)
df$Contamination <- as.numeric(df$Contamination)
df$Genome_Size <- as.numeric(df$Genome_Size)
df <- df[complete.cases(df), ]

# ==================== 2. 添加质量分组 ====================
df <- df %>%
  mutate(
    Quality_Group = case_when(
      Completeness > 90 & Contamination < 5  ~ "High-quality",
      Completeness >= 50 & Contamination < 10 ~ "Medium-quality",
      TRUE ~ "Low-quality"
    )
  )
df$Quality_Group <- factor(df$Quality_Group, 
                           levels = c("High-quality", "Medium-quality", "Low-quality"))

# ==================== 3. 定义颜色 ====================
quality_colors <- c("High-quality" = "#70C0AD",    
                    "Medium-quality" = "#EE9C75",  
                    "Low-quality" = "#9FABCC")

# ==================== 4. 绘制散点图 ====================
p1 <- ggplot(df, aes(x = Completeness, y = Contamination, 
                     color = Quality_Group, size = Genome_Size)) +
  geom_point(alpha = 0.8) +
  geom_smooth(size = 1.3, color = "#4c6389", method = "loess", se = TRUE) +   # 显示置信区间
  scale_color_manual(values = quality_colors) +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 1, colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.background = element_rect(colour = "black", size = 0.8)
  ) +
  labs(y = "Contamination (%)", x = "Completeness (%)", 
       color = "Quality", size = "Genome Size (bp)")

# 提取图例
p1_legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

# ==================== 5. 绘制边缘直方图（柱形） ====================
# 5.1 Completeness 直方图（顶部）
density_x <- ggplot(df, aes(x = Completeness, fill = Quality_Group)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", alpha = 1, bins = 30, color = "#969696") +
  scale_fill_manual(values = quality_colors) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank()
  ) +
  labs(y = "Density")

# 5.2 Contamination 直方图（右侧，翻转）
density_y <- ggplot(df, aes(x = Contamination, fill = Quality_Group)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", alpha = 1, bins = 30, color = "#969696") +
  scale_fill_manual(values = quality_colors) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_blank()
  ) +
  labs(x = "Density") +
  coord_flip()

# ==================== 6. 组合图形 ====================
combined <- ggdraw(
  plot_grid(
    plot_grid(density_x, p1, ncol = 1, rel_heights = c(0.3, 1)),
    plot_grid(p1_legend, density_y, ncol = 1, rel_heights = c(0.3, 1)),
    ncol = 2, rel_widths = c(1, 0.3)
  )
)

# 显示图形
combined
geom_density
# ==================== 7. 保存图片 ====================
ggsave("MAGs_completeness_contamination_quality_hist.pdf", combined, 
       width = 10, height = 8, units = "in", dpi = 300)