library(tidyverse)

# 1. 读取数据
rtab <- read.table("gene_presence_absence.Rtab", header=TRUE, row.names=1, check.names=FALSE)

n_genomes <- ncol(rtab)

# 2. 计算频率并定义基因分类
gene_freq <- data.frame(freq = (rowSums(rtab) / n_genomes) * 100) %>%
  mutate(category = case_when(
    freq >= 99 ~ "Core genes",
    freq >= 95 ~ "Soft core genes",
    freq >= 15 ~ "Shell genes",
    TRUE       ~ "Cloud genes"
  )) %>%
  # 将分类转为因子，确保图例顺序符合逻辑（从核心到外周）
  mutate(category = factor(category, levels = c("Cloud genes", "Shell genes", "Soft core genes", "Core genes")))

# 3. 绘图 A：不同分类不同颜色
plot_a <- ggplot(gene_freq, aes(x = freq, fill = category)) +
  # 使用 position = "stack" 确保分类颜色在直方图中堆叠显示
  geom_histogram(binwidth = 5, color = "white", boundary = 0, position = "stack") +
  # 自定义符合图片风格的颜色
  scale_fill_manual(values = c(
    "Core genes" = "#336456",      # 深蓝色
    "Soft core genes" = "#6699CC", # 浅蓝色
    "Shell genes" = "#99C",     # 浅绿色
    "Cloud genes" = "#339645"      # 深绿色
  )) +
  labs(
    x = "Gene frequency (%)", 
    y = "Number of genes", 
    title = "A",
    fill = "Gene Category"
  ) +
  theme_classic() +
  theme(legend.position = "right")

print(plot_a)