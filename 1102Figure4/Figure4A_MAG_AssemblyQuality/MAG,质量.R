# ==================== 加载包 ====================
library(dplyr)
library(ggplot2)
library(readr)

# ==================== 1. 读取丰度数据 ====================
abundance_file <- "abundance_mag.tsv" 
abundance <- read.delim(abundance_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(abundance)[1] <- "Genome"

abundance_mean <- abundance %>%
  mutate(Mean_Abundance = rowMeans(across(-Genome), na.rm = TRUE)) %>%
  select(Genome, Mean_Abundance)

# ==================== 2. 读取质量数据并定义分组 ====================
input_file <- "quality_report.tsv"
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(df)[1] <- "Genome" 
df <- df[, c("Genome", "Completeness", "Contamination", "Genome_Size")]

df$Completeness <- as.numeric(df$Completeness)
df$Contamination <- as.numeric(df$Contamination)
df$Genome_Size <- as.numeric(df$Genome_Size)
df <- df[complete.cases(df), ]

df <- df %>%
  mutate(
    Quality_Group = case_when(
      Completeness > 90 & Contamination < 5  ~ "High-quality",
      Completeness >= 50 & Contamination < 10 ~ "Medium-quality",
      TRUE ~ "Low-quality"
    )
  )
df$Quality_Group <- factor(df$Quality_Group, levels = c("High-quality", "Medium-quality", "Low-quality"))

# ==================== 3. 合并数据并过滤 ====================
plot_data <- merge(abundance_mean, df, by = "Genome")
plot_data <- plot_data %>% filter(Quality_Group != "Low-quality")
# 在合并数据后，画图前，加入这行代码：
plot_data <- plot_data %>% 
  filter(Mean_Abundance > 0) %>%
  filter(Quality_Group != "Low-quality") # 你原来的过滤保留
# ==================== 4. 计算 Spearman 相关性和文本标签 ====================
cor_test <- cor.test(plot_data$Mean_Abundance, plot_data$Completeness, method = "spearman")
R_val <- round(cor_test$estimate, 2)
P_val <- cor_test$p.value

# 处理 P 值，防止显示为 0，并且避免出现 "P=<" 的错误格式
if (P_val < 0.001) {
  label_text <- paste0("R=", R_val, ", P < 0.001")
} else {
  label_text <- paste0("R=", R_val, ", P=", sprintf("%.2f", P_val))
}

# ==================== 5. 定义颜色 ====================
quality_colors <- c("High-quality" = "#70C0AD",    
                    "Medium-quality" = "#EE9C75",  
                    "Low-quality" = "#9FABCC")

# ==================== 6. 绘制散点图 ====================
p <- ggplot(plot_data, aes(x = Mean_Abundance, y = Completeness, color = Quality_Group)) +
  geom_point(alpha = 0.6, size = 3) +  # 点的大小为3
  geom_smooth(method = "lm", color = "red", fill = "grey80", se = TRUE, linewidth = 1) +
  scale_color_manual(values = quality_colors) +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2),
                labels = parse(text = c("10^-3", "10^-2", "10^-1", "10^0", "10^1", "10^2"))) +
  scale_y_continuous(breaks = seq(60, 100, by = 20), limits = c(50, 105)) +
  
  # 将 R 和 P 值直接作为副标题（副标题默认在画框外部）
  labs(
    x = "Relative Abundance", 
    y = "Completeness (%)", 
    color = "MAG Quality",
    subtitle = label_text
  ) +
  
  # 使用 theme_bw() 获得完整的四条边框
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    
    # 【修改这里】副标题设置：靠左对齐，字体调大到 size = 8
    plot.subtitle = element_text(hjust = 0, size = 8, face = "plain"),
    plot.margin = margin(t = 25, r = 10, b = 10, l = 10) # 顶部留出足够空间
  )

# ==================== 7. 显示和保存 ====================
print(p)
ggsave("Figure_B_Complete.pdf", p, width = 6, height = 6)