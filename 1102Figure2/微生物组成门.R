library(tidyverse)
library(patchwork)
library(scales)

# ================= 1. Read abundance table =================
species <- read.delim(
  "bracken.P.0.2",
  row.names = 1,
  check.names = FALSE
)

# ================= 2. Read metadata =================
metadata <- read.delim(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# metadata format example:
# SampleID    Group
# S1          WM
# S2          NM
# S3          BM

group_df <- metadata %>%
  select(Samples, Group) %>%
  rename(Sample = Samples)

# Keep only samples existing in abundance table
group_df <- group_df %>%
  filter(Sample %in% colnames(species))

# Set group order
group_df$Group <- factor(
  group_df$Group,
  levels = c("WM", "NM", "BM")
)

# ================= 3. Select Top10 taxa =================
phylum_clean <- species %>%
  rownames_to_column("phylum") %>%
  mutate(total = rowSums(across(-phylum))) %>%
  arrange(desc(total)) %>%
  slice(1:10) %>%
  select(-total)

# ================= 4. Merge low abundance taxa =================
others <- species %>%
  filter(!rownames(.) %in% phylum_clean$phylum) %>%
  summarise(across(everything(), sum)) %>%
  mutate(phylum = "Others")

# ================= 5. Calculate relative abundance =================
df_phylum <- bind_rows(phylum_clean, others) %>%
  mutate(across(-phylum, ~ .x / sum(.x))) %>%
  pivot_longer(
    -phylum,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(group_df, by = "Sample")

# ================= 6. Sample order =================
target_phylum <- "Actinomycetota"

sample_order_df <- df_phylum %>%
  filter(phylum == target_phylum) %>%
  arrange(Group, desc(Abundance)) %>%
  pull(Sample)

df_phylum$Sample <- factor(
  df_phylum$Sample,
  levels = sample_order_df
)

# ================= 7. Taxa colors =================
color_map <- c(
  "Actinomycetota"  = "#3288bd",
  "Pseudomonadota"  = "#66c2a5",
  "Bacillota"       = "#9e0142",
  "Bacteroidota"    = "#fdae61",
  "Streptophyta"    = "#d53e4f",
  "Fusobacteriota"  = "#abdda4",
  "Euryarchaeota"   = "#f46d43",
  "Cyanobacteriota" = "#E9792F",
  "Apicomplexa"     = "#ffffbf",
  "Mycoplasmatota"  = "#F6C8A8",
  "Others"          = "#BDBDBD"
)

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
    Group = factor(
      group_name,
      levels = names(color_map_treat)
    )
  )
  
  ggplot(bar_data,
         aes(x = x,
             y = y,
             fill = Group)) +
    
    geom_tile() +
    
    scale_fill_manual(
      values = color_map_treat
    ) +
    
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        t = 5,
        r = 0,
        b = 5,
        l = 0
      ),
      legend.position = "none"
    )
}

# ================= 10. Area plot =================
plot_musk_area <- function(group_name,
                           data,
                           color_vec,
                           position = "middle") {
  
  sub_data <- data %>%
    filter(Group == group_name)
  
  sub_data$phylum <- factor(
    sub_data$phylum,
    levels = names(color_vec)
  )
  
  sub_data$x_numeric <- as.numeric(sub_data$Sample)
  
  p <- ggplot(
    sub_data,
    aes(
      x = x_numeric,
      y = Abundance,
      fill = phylum,
      group = phylum
    )
  ) +
    
    geom_area(
      position = "fill",
      alpha = 0.85,
      colour = NA
    ) +
    
    scale_fill_manual(
      values = color_vec,
      name = "Phylum"
    ) +
    
    scale_y_continuous(
      expand = c(0, 0),
      labels = percent
    ) +
    
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    
    theme_bw() +
    
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black"),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    ) +
    
    labs(
      x = NULL,
      y = NULL
    )
  
  if(position == "left") {
    
    p <- p +
      labs(y = "Relative Abundance") +
      
      theme(
        axis.line.y = element_line(
          colour = "black"
        )
      )
    
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
t_WM <- plot_top_bar(
  "WM",
  Treatments_colors
)

a_WM <- plot_musk_area(
  "WM",
  df_phylum,
  color_map,
  position = "left"
)

t_NM <- plot_top_bar(
  "NM",
  Treatments_colors
)

a_NM <- plot_musk_area(
  "NM",
  df_phylum,
  color_map,
  position = "middle"
)

t_BM <- plot_top_bar(
  "BM",
  Treatments_colors
)

a_BM <- plot_musk_area(
  "BM",
  df_phylum,
  color_map,
  position = "middle"
)

# ================= 12. Combine plots =================
col_WM <- t_WM / a_WM +
  plot_layout(
    heights = c(0.4, 20)
  )

col_NM <- t_NM / a_NM +
  plot_layout(
    heights = c(0.4, 20)
  )

col_BM <- t_BM / a_BM +
  plot_layout(
    heights = c(0.4, 20)
  )

final_plot <- (col_WM | col_NM | col_BM) +
  
  plot_layout(
    guides = "collect",
    widths = c(10, 9, 11)
  ) &
  
  theme(
    legend.position = "right",
    panel.spacing = unit(0, "cm")
  )

# ================= 13. Show plot =================
print(final_plot)

# ================= 14. Export PDF =================
ggsave(
  "Stacked_Bar_NoGap_Alpha0.8.pdf",
  final_plot,
  width = 12,
  height = 6,
  device = grDevices::cairo_pdf
)


##
# ================= 15. Export group mean relative abundance =================

library(dplyr)
library(tidyr)

# 每组平均相对丰度（基于样本平均）
group_mean <- df_phylum %>%
  group_by(Group, phylum) %>%
  summarise(
    Mean = mean(Abundance),
    SD = sd(Abundance),
    .groups = "drop"
  ) %>%
  mutate(
    Mean = Mean * 100,
    SD = SD * 100
  )

# 导出长表
write.csv(
  group_mean,
  "Group_mean_relative_abundance_long.csv",
  row.names = FALSE,
  quote = FALSE
)

# 导出宽表（论文推荐）
group_mean_wide <- group_mean %>%
  select(phylum, Group, Mean) %>%
  pivot_wider(
    names_from = Group,
    values_from = Mean
  )

write.csv(
  group_mean_wide,
  "Group_mean_relative_abundance_wide.csv",
  row.names = FALSE,
  quote = FALSE
)

# 查看结果
print(group_mean_wide)

# ================= 16. Overall relative abundance =================

overall_mean <- df_phylum %>%
  group_by(phylum) %>%
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

write.csv(
  overall_mean,
  "Overall_relative_abundance.csv",
  row.names = FALSE,
  quote = FALSE
)

print(overall_mean)