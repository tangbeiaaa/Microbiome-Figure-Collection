library(tidyverse)
library(scales)

# 1. 读取 Panaroo Rtab 文件
rtab <- read.table("gene_presence_absence.Rtab", header = TRUE, row.names = 1, check.names = FALSE)
pa_matrix <- as.matrix(rtab)
n_genomes <- ncol(pa_matrix)

# 2. 模拟随机抽样函数 (进行 100 次置换以产生平滑阴影)
n_iter <- 100
results <- list()

for (i in 1:n_iter) {
  shuffled <- pa_matrix[, sample(n_genomes)]
  pan_counts <- numeric(n_genomes)
  core_counts <- numeric(n_genomes)
  
  for (j in 1:n_genomes) {
    sub_mat <- shuffled[, 1:j, drop = FALSE]
    pan_counts[j] <- sum(rowSums(sub_mat) >= 1)     # 泛基因：并集
    core_counts[j] <- sum(rowSums(sub_mat) == j)    # 核心基因：交集
  }
  results[[i]] <- data.frame(Number = 1:n_genomes, Pan = pan_counts, Core = core_counts, Iter = i)
}

plot_data <- bind_rows(results) %>%
  pivot_longer(cols = c(Pan, Core), names_to = "Type", values_to = "Count")

# 3. 绘图
ggplot(plot_data, aes(x = Number, y = Count, color = Type, fill = Type, group = interaction(Type, Iter))) +
  # 绘制 100 条透明背景线形成阴影 (代表置信区间)
  geom_line(alpha = 0.05, size = 0.5) +
  # 绘制平均值趋势线 (Loess 平滑处理)
  stat_summary(aes(group = Type), fun = mean, geom = "line", size = 1.5) +
  # 设置配色
  scale_color_manual(values = c("Core" = "#C0392B", "Pan" = "#2471A3"),
                     labels = c("Core genome", "Pan genome")) +
  scale_fill_manual(values = c("Core" = "#C0392B", "Pan" = "#2471A3")) +
  # 坐标轴优化
  scale_x_continuous(breaks = pretty_breaks(10), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = comma) +
  labs(x = "Number of Genomes", 
       y = "Number of Gene Families",
       title = "Pangenome and Core genome Rarefaction Curves") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")) +
  # 修正图例透明度
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))