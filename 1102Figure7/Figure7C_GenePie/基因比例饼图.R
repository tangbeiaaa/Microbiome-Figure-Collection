# 加载必要的库
library(ggplot2)
library(dplyr)

# 1. 准备数据
data <- data.frame(
  category = c("Core genes", "Soft core genes", "Shell genes", "Cloud genes"),
  count = c(1599, 98, 225, 219)
)

# 计算百分比
data <- data %>%
  mutate(percentage = count / sum(count) * 100,
         label = paste0(category, "\n(", count, ", ", round(percentage, 2), "%)"))

# 2. 绘图
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # 去除背景坐标轴
  scale_fill_brewer(palette = "Set3") + # 设置配色方案
  labs(title = "Pangenome Composition",
       fill = "Gene Categories") +
  geom_text(aes(label = ifelse(percentage > 5, label, "")), # 仅在大区域显示文本避免重叠
            position = position_stack(vjust = 0.5))