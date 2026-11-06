library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(gridExtra)
library(patchwork)
library(readxl)
# 读取数据
quality_data <- read_excel("quality_report_MAGs.xlsx")

# 查看数据结构
head(quality_data)
str(quality_data)

# 使用严格标准定义质量等级
quality_data <- quality_data %>%
  mutate(
    Quality_Group = case_when(
      # 高质量: Completeness > 90% 且 Contamination < 5%
      Completeness > 90 & Contamination < 5 ~ "High-quality",
      # 中等质量: Completeness >= 50% 且 Contamination < 10%
      Completeness >= 50 & Contamination < 10 ~ "Medium-quality",
      # 低质量: 其他情况
      TRUE ~ "Low-quality"
    )
  )

# 按 Source 分组计算质量分布
quality_summary_by_source <- quality_data %>%
  group_by(Source, Quality_Group) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Source) %>%
  mutate(
    Total = sum(Count),
    Percentage = Count / Total * 100,
    Label = paste0(Quality_Group, "\n", round(Percentage, 1), "%")
  ) %>%
  ungroup()

# 创建颜色映射
quality_colors <- c(
  "High-quality" = "#70C0AD",    
  "Medium-quality" = "#EE9C75",  
  "Low-quality" = "#9FABCC"      
)

# 创建饼状图函数，图例在底部
create_pie_chart <- function(source_data, source_name) {
  ggplot(source_data, aes(x = "", y = Count, fill = Quality_Group)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = quality_colors, name = "") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),  # 小方块
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = paste(source_name, "MAGs"),
      subtitle = paste("Total:", unique(source_data$Total), "MAGs")
    ) +
    geom_text(
      aes(label = Label), 
      position = position_stack(vjust = 0.5),
      size = 3.5, 
      color = "white", 
      fontface = "bold"
    )
}

# 获取每个来源的数据
mmbc_data <- quality_summary_by_source %>% filter(Source == "Study samples")
published_data <- quality_summary_by_source %>% filter(Source == "Public datasets")

# 创建饼状图
p_mmbc <- create_pie_chart(mmbc_data, "Study samples")
p_published <- create_pie_chart(published_data, "Public datasets")

# 使用patchwork组合图形，并收集图例到底部
combined_plot <- (p_mmbc | p_published) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')

# 显示最终图形
print(combined_plot)

# 保存图形
ggsave("genome_quality_by_source.png", combined_plot, width = 10, height = 6, dpi = 300)