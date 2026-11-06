# 加载必备包
library(ggplot2)
library(dplyr)
library(cowplot)

# ------------------ 1. 读取数据 ------------------
abund <- read.table("abundance.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tax_anno <- read.table("tax.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 你的分组信息（保持原样）
wm_cols <- c("A10_1","A11_1","A12_1","A7_1","A8_1","A9_1","B7_1","C10_1","C11_1","C12_1","C7_1","C8_1","C9_1")
nm_cols <- c("A13_1","A15_1","A16_1","A18_1","B13_1","B14_1","B15_1","B17_1","B18_1","C13_1")
bm_cols <- c("A1_1","A3_1","A4_1","A5_1","A6_1","B1_1","B2_1","B3_1","B4_1","B5_1","B6_1","C1_1","C2_1","C3_1","C4_1","C5_1","C6_1")
all_sample_cols <- c(wm_cols, nm_cols, bm_cols)

# ================= 核心修改：SGB 聚合为 属 =================

# 1. 建立 SGB_ID 到 属名 的映射字典，并去除 "g__" 前缀
id_to_genus <- setNames(tax_anno$Genus, tax_anno$ID)
clean_genus <- gsub("^g__", "", id_to_genus)               # 去掉前缀
# 🔴 修改1：把空字符串替换为你要求的 "unmapped"
clean_genus[clean_genus == ""] <- "unmapped"           

# 2. 为 abund 的每一行匹配对应的属名
matched_genus <- clean_genus[rownames(abund)]
matched_genus[is.na(matched_genus)] <- rownames(abund)[is.na(matched_genus)]

# 3. 把所有 SGB 丰度按照属名汇总 (聚合为属水平丰度)
genus_abund <- rowsum(abund[, all_sample_cols, drop = FALSE], group = matched_genus)

# 4. 计算所有属的平均相对丰度，找出 Top20 属
row_means <- rowMeans(genus_abund[, all_sample_cols], na.rm = TRUE)
top20_genus_names <- names(sort(row_means, decreasing = TRUE))[1:20]
top20_abund <- genus_abund[top20_genus_names, ]

# ================= 至此，top20_abund 就是丰度前20的属了 =================

# ------------------ 3. 计算各组流行率 ------------------
calc_prev <- function(data, group_cols) {
  present_counts <- rowSums(data[, group_cols, drop = FALSE] > 0)
  return(present_counts / length(group_cols))
}

wm_prev <- calc_prev(top20_abund, wm_cols)
nm_prev <- calc_prev(top20_abund, nm_cols)
bm_prev <- calc_prev(top20_abund, bm_cols)

# ------------------ 4. 整理画图数据（Y轴标签取属名） ------------------
plot_data <- data.frame(
  Genus = rep(rownames(top20_abund), 3),
  Group = c(rep("WM", 20), rep("NM", 20), rep("BM", 20)),
  Prevalence = c(wm_prev, nm_prev, bm_prev)
)

# 🔴 修改2：自定义 Factor 顺序，将 "unmapped" 强制放在最下面
top20_genera <- rownames(top20_abund)
if ("unmapped" %in% top20_genera) {
  # 如果 unmapped 在 Top20 中，将其拿出来放在 levels 的【第1位】（即图的最底部）
  top20_rest <- top20_genera[top20_genera != "unmapped"]
  plot_data$Genus <- factor(plot_data$Genus, 
                            levels = c("unmapped", rev(top20_rest)))
} else {
  # 如果 Top20 里没有 unmapped，直接按丰度从高到低排序（即图的最底部是最低丰度属）
  plot_data$Genus <- factor(plot_data$Genus, 
                            levels = rev(top20_genera))
}

# ------------------ 5. 绘制三个分组子图 ------------------
groups <- unique(plot_data$Group)
colors <- c("WM"="#dfc27d","NM" = "#bf812d","BM" = "#543005")

plots <- lapply(seq_along(groups), function(i) {
  g <- groups[i]
  dat <- subset(plot_data, Group == g)
  
  p <- ggplot(dat, aes(x = Genus, y = Prevalence)) +
    geom_bar(stat = "identity", fill = colors[g]) +
    coord_flip() +
    labs(title = g) +  
    theme_classic() +
    theme(
      text = element_text(family = 'sans', size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold") 
    ) +
    scale_y_continuous(expand = expansion(0, 0), limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.50, 0.75, 1.0)) +
    geom_hline(yintercept = c(0.25, 0.50, 0.75, 1.0), colour = 'grey60', lwd = 0.3, linetype = "dashed")
  
  # 只保留第1个图（WM）的 Y 轴物种标签
  if (i == 1) {
    p <- p + theme(axis.text.y = element_text(size = 9, color = "black", hjust = 1))
  }
  return(p)
})

# ------------------ 6. 拼图并导出 ------------------
final_plot <- plot_grid(plotlist = plots, ncol = 3, rel_widths = c(1.5, 1, 1), align = "h")
print(final_plot)

if(!dir.exists("results")) dir.create("results")
pdf("results/Top20_Genera_Prevalence.pdf", height = 7.2, width = 8)
final_plot
dev.off()