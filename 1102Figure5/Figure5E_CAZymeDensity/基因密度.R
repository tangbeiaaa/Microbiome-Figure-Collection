#####
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
#kegg_gene.count
#gene.count
# ==========================================
# 1. 读取原始数据
# ==========================================
df_anno <- read.delim("annotation2self.txt", header = TRUE, stringsAsFactors = FALSE)
df_caz <- read.table("kegg_gene.count", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

# ==========================================
# 2. 提取丰度列，并保留这些列用于后续聚合
# ==========================================
# 自动识别样本列（列名以 字母+数字+下划线 开头，如 A10_1、B13_1 等）
sample_cols <- names(df_anno)[grep("^[A-Z][0-9]+_", names(df_anno))]

# 提取核心列（**必须把样本列 all_of(sample_cols) 带进来**）
df_info <- df_anno %>%
  select(
    ID, 
    Genus, 
    Genome_Size = size, 
    Completeness = completeness, 
    Contamination = contamination,
    all_of(sample_cols)  # 🔴 保留原始样本列，便于后续按属聚合
  )

# ==========================================
# 3. 合并数据：将 CAZyme 数量追加到主表
# ==========================================
df_merge <- df_info %>%
  inner_join(df_caz, by = "ID")

# ==========================================
# 4. 数据过滤与清洗（过滤低质量，排除 unmapped g__）
# ==========================================
df_clean <- df_merge %>%
  filter(Genus != "g__") %>%
  filter(Completeness > 50 & Contamination < 10)

# ==========================================
# 5. 核心计算：CAZyme 密度（保留每个 MAG 的密度）
# ==========================================
df_clean <- df_clean %>%
  mutate(CAZyme_density = KEGG / (Genome_Size / 1000000))
#df_clean <- df_clean %>%
#  mutate(CAZyme_density = Cazy / (Genome_Size / 1000000))

# ==========================================
# 6. 筛选 Top 20 属：【改为按群落总相对丰度（sum）排序】
# ==========================================
# 6.1 在属水平上，将同属所有 MAG 的绝对丰度相加
genus_aggregated <- df_clean %>%
  group_by(Genus) %>%
  summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE))) %>%
  mutate(
    # 6.2 计算该属在所有样本中的平均总相对丰度
    Total_Abundance = rowMeans(across(all_of(sample_cols)), na.rm = TRUE)
  ) %>%
  arrange(desc(Total_Abundance))

# 6.3 提取总丰度最高的前 20 个属名
top20_genus <- genus_aggregated %>%
  slice_head(n = 20) %>%
  pull(Genus)

# ==========================================
# 7. 准备绘图数据 & 锁定 X 轴顺序
# ==========================================
df_plot <- df_clean %>%
  filter(Genus %in% top20_genus) %>%
  mutate(
    # 去掉 g__ 前缀，让图上的标签更干净
    Genus = gsub("g__", "", Genus) 
  ) %>%
  # 确保排序与 top20_genus 严格一致 (注意 top20_genus 也要去掉前缀)
  mutate(Genus = factor(Genus, levels = gsub("g__", "", top20_genus)))

# ==========================================
# 8. 绘制箱线图
# ==========================================
# 生成 20 种柔和配色
colors_20 <- grDevices::hcl.colors(length(unique(df_plot$Genus)), palette = "Set 3")

p <- ggplot(df_plot, aes(x = Genus, y = CAZyme_density, fill = Genus)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) + 
  theme_bw() + 
  labs(
    y = "CAZymes per 1 M genome",
    x = NULL, 
    title = "CAZyme density of the top 20 genera (by total abundance)" # 修改了标题，呼应排序逻辑
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 10, face = "italic"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none", 
    panel.grid.major.x = element_blank() 
  ) +
  scale_fill_manual(values = colors_20)

# ==========================================
# 9. 查看并保存图片
# ==========================================
print(p)

