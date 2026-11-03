# =========================
# 1. 加载包
# =========================
library(dplyr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(readxl)

# =========================
# 2. 读取数据
# =========================
df <- read_excel("MAX_ANI_G.xlsx")

# 处理非数值
df$MAX_ani[df$MAX_ani == "N/A"] <- NA
df$MAX_ani <- as.numeric(df$MAX_ani)

df <- df %>%
  filter(!is.na(MAX_ani))

# =========================
# 3. 选择前 N 个 genus
# =========================
N <- 10

genus_counts <- df %>%
  count(genus, sort = TRUE)

top_genera <- genus_counts$genus[1:N]

df_top <- df %>%
  filter(genus %in% top_genera)

# =========================
# 4. 定义 Known / Novel
# =========================
df_top <- df_top %>%
  mutate(
    Status = case_when(
      MAX_ani >= 95 ~ "Known",
      MAX_ani < 83 ~ "Novel",
      TRUE ~ "Intermediate"
    )
  )

# =========================
# 5. genus 排序（按 median ANI）
# =========================
genus_order <- df_top %>%
  group_by(genus) %>%
  summarise(median_ani = median(MAX_ani)) %>%
  arrange(median_ani) %>%
  pull(genus)

df_top$genus <- factor(df_top$genus,
                       levels = genus_order)

# =========================
# 6. 左侧 ridge plot
# =========================
ridge_plot <- ggplot(
  df_top,
  aes(
    x = MAX_ani,
    y = genus,
    fill = after_stat(x)
  )
) +
  
  geom_density_ridges_gradient(
    scale = 1.4,
    rel_min_height = 0.01,
    alpha = 0.9
  ) +
  
  scale_fill_gradientn(
    colors = c("#2E7AB8", "white", "#80D3A9")
  ) +
  
  # Known / Novel 阈值线
  geom_vline(
    xintercept = 95,
    linetype = "dashed",
    color = "#80D3A9",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = 83,
    linetype = "dashed",
    color = "#2E7AB8",
    linewidth = 0.8
  ) +
  
  # 标注
  annotate(
    "text",
    x = 96,
    y = length(genus_order) + 0.8,
    label = "Known",
    color = "#80D3A9",
    size = 4
  ) +
  
  annotate(
    "text",
    x = 84,
    y = length(genus_order) + 0.8,
    label = "Novel",
    color = "#2E7AB8",
    size = 4
  ) +
  
  scale_x_continuous(
    limits = c(70, 100),
    breaks = seq(70, 100, 5)
  ) +
  
  labs(
    x = "MAX ANI (%)",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.position = "none"
  )

# =========================
# 7. 右侧双向柱状图
# =========================

bar_data <- df_top %>%
  filter(Status != "Intermediate") %>%
  group_by(genus, Status) %>%
  summarise(Count = n(), .groups = "drop")

# Novel 向左
bar_data <- bar_data %>%
  mutate(
    Count = ifelse(Status == "Novel",
                   -Count,
                   Count)
  )

# 绘图
bar_plot <- ggplot(
  bar_data,
  aes(
    x = Count,
    y = genus,
    fill = Status
  )
) +
  
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      "Known" = "#80D3A9",
      "Novel" = "#2E7AB8"
    )
  ) +
  
  scale_x_continuous(
    labels = abs,
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Number of SGBs",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.position = "top"
  )

# =========================
# 8. 拼图
# =========================
combined_plot <- ridge_plot + bar_plot +
  plot_layout(widths = c(2.2, 1))

# 显示
print(combined_plot)

# 保存
# ggsave("ANI_known_novel.pdf",
#        combined_plot,
#        width = 10,
#        height = 6,
#        dpi = 300)