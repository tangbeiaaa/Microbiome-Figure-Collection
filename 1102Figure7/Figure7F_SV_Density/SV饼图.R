# 加载必要的库
library(ggplot2)
library(dplyr)

# 1. 准备数据
data <- data.frame(
  category = c("transloc", "inversion", "deletion", "insertion", "duplication"),
  count = c(4, 4, 166, 85, 6)
)

# 计算百分比
data <- data %>%
  mutate(percentage = count / sum(count) * 100,
         label = paste0(category, "\n(", count, ", ", round(percentage, 2), "%)"))

# 自定义颜色（顺序不影响，会自动按category匹配）
type_colors <- c("insertion" = "#E2C56A",
                 "deletion" = "#A35338",
                 "transloc" = "#BAC8DC",
                 "inversion" = "#EAB7BF",
                 "duplication" = "#57837C",
                 "contraction"= "#9a9a9a")

# 绘制圆环图（使用自定义颜色）
ggplot(data, aes(x = 2, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  theme_void() +
  scale_fill_manual(values = type_colors) +   # 替换为自定义颜色
  labs(title = "SV Composition", fill = "type") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5))


###
library(ggplot2)
library(dplyr)

# 读取数据（假定文件已统一为制表符分隔）
data <- read.table("avg_density_by_group.tsv", header = TRUE, sep = "\t")

# 重命名列
colnames(data) <- c("Group", "Type", "Mean", "SD")  # SD 列仍保留但不用

data <- data %>%
  complete(Type, Group,
           fill = list(Mean = 0,
                       SD = 0))

# -----------------------------
# 重新设置 factor（必须在complete后）
# -----------------------------
data$Group <- factor(data$Group,
                     levels = c("WM", "NM", "BM"))

# 按 SV 类型的总体平均密度从高到低排序横坐标
type_order <- data %>%
  group_by(Type) %>%
  summarise(OverallMean = mean(Mean, na.rm = TRUE)) %>%
  arrange(desc(OverallMean)) %>%
  pull(Type)

data$Type <- factor(data$Type, levels = type_order)
# 自定义颜色
my_colors <- c("WM" = "#f6e8c3", "NM" = "#bf812d", "BM" = "#543005")

# 绘图（无误差棒）
p <- ggplot(data, aes(x = Type, y = Mean, fill = Group)) +
  
  geom_bar(stat = "identity",
           position = position_dodge2(width = 0.9,
                                      preserve = "single"),
           na.rm = TRUE) +
  
  labs(x = "SV type",
       y = "Mean SV density (per Mb)",
       title = "Structural variant density by group") +
  
  scale_fill_manual(values = my_colors) +
  
  theme_minimal() +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

p<- p +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(p)

######堆叠著柱状图
# 加载必要的包
library(ggplot2)
library(dplyr)

data <- read.table("sv_group_type_counts.txt", header = FALSE, sep = "\t", 
                   col.names = c("Group", "Type", "Count"))

data$Type <- as.character(data$Type)

# 设置横坐标顺序为 WM, NM, BM
data$Group <- factor(data$Group, levels = c("WM", "NM", "BM"))

data <- data %>%
  group_by(Group) %>%
  mutate(Total = sum(Count),
         Percent = Count / Total * 100)

type_colors <- c("insertion" = "#BEBADA",
                 "deletion" = "#8DD3C7",
                 "transloc" = "#80B1D3",
                 "inversion" = "#FB8072",
                 "duplication" = "#FFFFB3",
                 "contraction"= "#9a9a9a")

present_types <- unique(data$Type)
color_palette <- type_colors[names(type_colors) %in% present_types]

p <- ggplot(data, aes(x = Type, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(x = "SV type", y = "Number of SV events",
       title = "Structural variant counts by group") +
  scale_fill_manual(values = c("WM" = "#f6e8c3", "NM" = "#bf812d", "BM" = "#543005")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

print(p)
