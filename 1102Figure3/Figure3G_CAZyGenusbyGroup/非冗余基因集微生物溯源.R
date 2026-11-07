################################################xin
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)
# ============================================================
# 1. 读取 all_gene.list，提取 ORF → CAZy（一级，不剔除编号0）
# ============================================================
orf_cazy <- fread("all_gene.list", header = TRUE, sep = "\t",
                  fill = TRUE, na.strings = "")
colnames(orf_cazy)[1] <- "ORF"   # Name → ORF

orf_cazy_long <- orf_cazy[, .(ORF, CAZy)] %>%
  as_tibble() %>%
  filter(!is.na(CAZy), CAZy != "", CAZy != "-") %>%
  mutate(CAZy = trimws(CAZy)) %>%
  separate_rows(CAZy, sep = ",") %>%      # 兼容逗号分隔多值
  mutate(CAZy = trimws(CAZy)) %>%
  filter(grepl("^(GH|GT|PL|CE|AA|CBM)\\d+", CAZy)) %>%
  # 一级分类：GH0/GH33 → GH，CBM48 → CBM；编号0条目保留
  mutate(CAZy = sub("\\d+.*$", "", CAZy)) %>%
  distinct(ORF, CAZy)

cat("ORF-CAZy(一级) 映射条数：", nrow(orf_cazy_long), "\n")
cat("一级 CAZy 类别：", paste(sort(unique(orf_cazy_long$CAZy)), collapse = ", "), "\n")

# ============================================================
# 2. ORF 丰度
# ============================================================
ORF_ab <- fread("all_gene.TPM", check.names = FALSE)
ORF_ab_long <- melt(ORF_ab, id.vars = "Name",
                    variable.name = "Samples",
                    value.name = "Abundance")
setnames(ORF_ab_long, "Name", "ORF")

# ============================================================
# 3. ORF → 分类
# ============================================================
orf_tax <- fread("all_nucleotide.tax", header = TRUE, sep = "\t",
                 fill = TRUE, na.strings = "")
setnames(orf_tax, "Name", "ORF")

# ============================================================
# 4. 合并
# ============================================================
# ============================================================
# 4. 合并（优化顺序，先缩小再 join）
# ============================================================

# 4.1 先把 orf_tax 自身过滤干净，只留必要列
orf_tax_clean <- orf_tax %>%
  filter(!is.na(taxid), taxid != "", taxid != "Unassigned", taxid != "-") %>%
  filter(!Genus %in% c("Homo", "Chordata", "Babesia", "Mus", "Rattus")) %>%
  select(ORF, taxid, Genus, Phylum, Class, Order, Family)  # 按需保留，列越少越省内存
gc()

# 4.2 ORF 丰度先过滤非零
ORF_ab_long <- ORF_ab_long %>% filter(Abundance > 0)
gc()

# 4.3 关键：先 inner_join 较小的 orf_cazy_long，立刻缩小行数
df <- ORF_ab_long %>%
  inner_join(orf_cazy_long, by = "ORF") %>%   # ← 放最前，行数骤降
  inner_join(orf_tax_clean, by = "ORF")        # ← 再挂 taxonomy
gc()

cat("合并后数据行数：", nrow(df), "\n")

# ============================================================
# 5. 每个一级 CAZy 的属水平贡献
# ============================================================
cazy_contrib <- df %>%
  group_by(CAZy, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(CAZy) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  arrange(CAZy, desc(Prop)) %>%
  ungroup()

fwrite(as.data.table(cazy_contrib),
       "CAZy_L1_genus_contribution_full.tsv",
       sep = "\t", quote = FALSE)

# 一级 CAZy 按总丰度排序（全部 6 类）
top_cazy <- cazy_contrib %>%
  group_by(CAZy) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(CAZy)

cat("一级 CAZy 按总丰度排序：\n")
print(top_cazy)

# ============================================================
# 6. 饼图：全局属贡献
# ============================================================
df_genus_total <- cazy_contrib %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  mutate(Genus = case_when(
    is.na(Genus) | Genus == "" ~ "Unassigned",
    Genus %in% c("Unassigned", "Unclassified") ~ "Unassigned",
    TRUE ~ Genus
  )) %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  arrange(desc(Prop))

top_n <- 10
top_genera <- df_genus_total %>%
  filter(Genus != "Unassigned") %>%
  slice_max(Prop, n = top_n) %>%
  pull(Genus)

base_colors <- c(
  "Corynebacterium" = "#3188ba",
  "Pseudomonas"     = "#65bfa4",
  "Pseudomonas_E"   = "#fbae61",
  "Bacteroides"     = "#f46d43",
  "Atopostipes"     = "#3D9F3C",
  "Mammaliicoccus"  = "#9d0e43",
  "Burkholderia"    = "#d53e4f",
  "Dietzia"         = "#abd7a2",
  "Enterococcus"    = "#A1CDEC",
  "Alloiococcus"    = "#f1c8a8",
  "Nitriliruptor"   = "#e9792e",
  "Streptococcus"   = "#f5f8be",
  "Shewanella"      = "#c7b3d6",
  "Staphylococcus"  = "#fddbc7",
  "Neisseria"       = "#c5d9e8",
  "Nostoc"          = "#b8e186",
  "Arthrobacter"    = "#f4a582",
  "Enterobacter"    = "#92c5de",
  "Providencia"     = "#d6604d",
  "Acinetobacter"   = "#4393c3",
  "Others"          = "#bdbdbd",
  "Unassigned"      = "#636363"
)

# ---- 图1：含 Others + Unassigned ----
df_pie_all <- df_genus_total %>%
  mutate(Group = case_when(
    Genus == "Unassigned" ~ "Unassigned",
    Genus %in% top_genera ~ Genus,
    TRUE ~ "Others"
  )) %>%
  group_by(Group) %>%
  summarise(Prop = sum(Prop), .groups = "drop") %>%
  arrange(desc(Prop))

group_order <- c(
  setdiff(df_pie_all$Group, c("Unassigned", "Others")),
  intersect(c("Others", "Unassigned"), df_pie_all$Group)
)
df_pie_all$Group <- factor(df_pie_all$Group, levels = group_order)
df_pie_all <- df_pie_all %>%
  mutate(label = paste0(Group, "\n", percent(Prop, accuracy = 0.1)))

extra_genera <- setdiff(levels(df_pie_all$Group), names(base_colors))
if (length(extra_genera) > 0) {
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(extra_genera))
  names(extra_colors) <- extra_genera
  base_colors <- c(base_colors, extra_colors)
}
pie_colors_all <- base_colors[levels(df_pie_all$Group)]

p_pie_all <- ggplot(df_pie_all, aes(x = "", y = Prop, fill = Group)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pie_colors_all) +
  theme_void() +
  labs(fill = "Genus",
       title = "CAZy 属贡献总览（一级，含 Others / Unassigned）")

print(p_pie_all)
ggsave("Pie_CAZy_L1_Genus_All.pdf", p_pie_all, width = 8, height = 8)
ggsave("Pie_CAZy_L1_Genus_All.png", p_pie_all, width = 8, height = 8, dpi = 300)

# ---- 图2：去掉 Others + Unassigned ----
df_pie_clean <- df_genus_total %>%
  filter(!Genus %in% c("Unassigned", "Others")) %>%
  slice_max(Prop, n = top_n) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  arrange(desc(Prop))
df_pie_clean$Genus <- factor(df_pie_clean$Genus, levels = df_pie_clean$Genus)
df_pie_clean <- df_pie_clean %>%
  mutate(label = paste0(Genus, "\n", percent(Prop, accuracy = 0.1)))

pie_colors_clean <- base_colors[levels(df_pie_clean$Genus)]
extra_genera2 <- setdiff(levels(df_pie_clean$Genus), names(pie_colors_clean))
if (length(extra_genera2) > 0) {
  extra_colors2 <- colorRampPalette(brewer.pal(12, "Set3"))(length(extra_genera2))
  names(extra_colors2) <- extra_genera2
  pie_colors_clean <- c(pie_colors_clean, extra_colors2)
}
pie_colors_clean <- pie_colors_clean[levels(df_pie_clean$Genus)]

p_pie_clean <- ggplot(df_pie_clean, aes(x = "", y = Prop, fill = Genus)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pie_colors_clean) +
  theme_void() +
  labs(fill = "Genus",
       title = "CAZy 属贡献（一级，去 Others / Unassigned）")

print(p_pie_clean)
ggsave("Pie_CAZy_L1_Genus_Clean.pdf", p_pie_clean, width = 8, height = 8)
ggsave("Pie_CAZy_L1_Genus_Clean.png", p_pie_clean, width = 8, height = 8, dpi = 300)

# ============================================================
# 7. 堆叠柱状图：一级 CAZy（GH/GT/PL/CE/AA/CBM）
# ============================================================
df_clean <- cazy_contrib %>%
  filter(!Genus %in% c("Others", "Unassigned", "Unclassified"),
         !is.na(Genus), Genus != "") %>%
  filter(CAZy %in% top_cazy)

top10_genera <- df_clean %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

df_plot <- df_clean %>% filter(Genus %in% top10_genera)

# 每个一级 CAZy 内重新归一化
df_plot <- df_plot %>%
  group_by(CAZy) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  ungroup()

# CAZy 排序：按 Corynebacterium 贡献降序；若不存在则按总贡献排序
target_genus <- "Corynebacterium"
if (target_genus %in% df_plot$Genus) {
  cazy_order <- df_plot %>%
    filter(Genus == target_genus) %>%
    group_by(CAZy) %>%
    summarise(target_contrib = sum(Prop), .groups = "drop") %>%
    arrange(desc(target_contrib)) %>%
    pull(CAZy)
} else {
  cazy_order <- df_plot %>%
    group_by(CAZy) %>%
    summarise(t = sum(Prop), .groups = "drop") %>%
    arrange(desc(t)) %>%
    pull(CAZy)
}
cazy_order <- c(cazy_order, setdiff(unique(df_plot$CAZy), cazy_order))
df_plot$CAZy <- factor(df_plot$CAZy, levels = cazy_order)

# 属排序：Corynebacterium 固定底部
genus_order <- df_plot %>%
  group_by(Genus) %>%
  summarise(total = sum(Prop), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Genus)
if (target_genus %in% genus_order) {
  genus_order <- c(target_genus, setdiff(genus_order, target_genus))
}
df_plot$Genus <- factor(df_plot$Genus, levels = genus_order)

extra_genera <- setdiff(levels(df_plot$Genus), names(base_colors))
if (length(extra_genera) > 0) {
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(extra_genera))
  names(extra_colors) <- extra_genera
  base_colors <- c(base_colors, extra_colors)
}
mycol <- base_colors[levels(df_plot$Genus)]

p_bar <- ggplot(df_plot, aes(x = CAZy, y = Prop, fill = Genus)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(
    x = "CAZy class (level 1)",
    y = "Relative contribution (%)",
    fill = "Genus",
    title = "一级 CAZy 微生物贡献堆叠图"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  )

print(p_bar)
##########################

#######################################fenzu

library(patchwork) # 确保加载了拼接包

# ============================================================
# 4.4 补充：挂载分组信息（新加）
# ============================================================
# 读取 metadata 并设定因子顺序
group_info <- read.delim("metadata.txt", header = TRUE, stringsAsFactors = FALSE)
group_info$Group <- factor(group_info$Group, levels = c("WM", "NM", "BM"))

# 将分组信息 join 到 df 上
df_grp <- df %>%
  left_join(group_info %>% select(Samples, Group), by = "Samples") %>%
  filter(!is.na(Group))

cat("挂载分组后数据行数：", nrow(df_grp), "\n")

# ============================================================
# 5. 每个一级 CAZy 的属水平贡献（按组）
# ============================================================
cazy_contrib_grp <- df_grp %>%
  group_by(Group, CAZy, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, CAZy) %>%
  filter(sum(Abundance) > 0) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  ungroup()

# 导出数据表
fwrite(as.data.table(cazy_contrib_grp),
       "CAZy_L1_genus_contribution_byGroup.tsv",
       sep = "\t", quote = FALSE)

# ============================================================
# 6. 取 Top 10 属（跨所有组合并，保证三组颜色/图例一致）
# ============================================================
df_clean <- cazy_contrib_grp %>%
  filter(!Genus %in% c("Others", "Unassigned", "Unclassified"),
         !is.na(Genus), Genus != "")

top10_genera <- df_clean %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

df_plot_grp <- df_clean %>% filter(Genus %in% top10_genera)

# 每个组内、每个 CAZy 重新归一化
df_plot_grp <- df_plot_grp %>%
  group_by(Group, CAZy) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  ungroup()

# ============================================================
# 7. 颜色方案（沿用你原来的配色）
# ============================================================
base_colors <- c(
  "Corynebacterium" = "#3188ba",
  "Pseudomonas_E"   = "#fbae61", "Bacteroides" = "#f46d43",
  "Atopostipes"     = "#3D9F3C", "Mammaliicoccus" = "#9d0e43",
  "Burkholderia"    = "#d53e4f", "Dietzia" = "#abd7a2",
  "Enterococcus"    = "#A1CDEC", "Alloiococcus" = "#f1c8a8",
  "Nitriliruptor"   = "#e9792e", "Streptococcus" = "#f5f8be",
  "Shewanella"      = "#c7b3d6", "Staphylococcus" = "#fddbc7",
  "Neisseria"       = "#c5d9e8", "Nostoc" = "#b8e186",
  "Arthrobacter"    = "#f4a582", "Enterobacter" = "#92c5de",
  "Providencia"     = "#d6604d", "Acinetobacter" = "#4393c3","Trueperella"     = "#66c2a5",
  "Pseudomonas"     = "#abdda4",
  "Others" = "#bdbdbd", "Unassigned" = "#636363"
)

# 补齐缺失属的颜色
extra_genera <- setdiff(unique(as.character(df_plot_grp$Genus)), names(base_colors))
if (length(extra_genera) > 0) {
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(extra_genera))
  names(extra_colors) <- extra_genera
  base_colors <- c(base_colors, extra_colors)
}
mycol_grp <- base_colors[unique(as.character(df_plot_grp$Genus))]
mycol_grp <- mycol_grp[!is.na(names(mycol_grp))]

# 设定排序目标
target_genus <- "Corynebacterium"

# ============================================================
# 8. 分组画堆叠图（每组内独立按 target_genus 的 Prop 降序）
# ============================================================
plot_list <- list()

for (g in levels(df_plot_grp$Group)) {
  
  d <- df_plot_grp %>% filter(Group == g)
  if (nrow(d) == 0) next
  
  # 计算当前组内，每条 CAZy 中 target_genus 的相对比例 (Prop)
  if (target_genus %in% unique(as.character(d$Genus))) {
    ord_df <- d %>%
      mutate(Genus_chr = as.character(Genus)) %>%
      filter(Genus_chr == target_genus) %>%
      group_by(CAZy) %>%
      summarise(prop_val = sum(Prop, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(prop_val))
    
    cazy_order <- c(ord_df$CAZy, setdiff(unique(as.character(d$CAZy)), ord_df$CAZy))
  } else {
    # 兜底：如果没有 Corynebacterium，按总丰度排
    cazy_order <- d %>%
      group_by(CAZy) %>%
      summarise(t = sum(Abundance), .groups = "drop") %>%
      arrange(desc(t)) %>% pull(CAZy)
  }
  
  # 设定当前组的因子水平
  d$CAZy <- factor(as.character(d$CAZy), levels = cazy_order)
  
  # 属的堆叠顺序：target_genus 固定底部，其余按丰度降序
  genus_ord <- d %>%
    group_by(Genus) %>%
    summarise(total = sum(Prop), .groups = "drop") %>%
    arrange(desc(total)) %>% pull(Genus) %>% as.character()
  if (target_genus %in% genus_ord) {
    genus_ord <- c(target_genus, setdiff(genus_ord, target_genus))
  }
  d$Genus <- factor(as.character(d$Genus), levels = genus_ord)
  
  # 画单张图
  p <- ggplot(d, aes(x = CAZy, y = Prop, fill = Genus)) +
    geom_col(width = 0.85, color = "white", linewidth = 0.15) +
    scale_fill_manual(values = mycol_grp) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(
      x = NULL,
      y = if(g == "WM") "Relative contribution (%)" else NULL,
      title = g
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  plot_list[[g]] <- p
}

# ============================================================
# 9. 拼图并保存
# ============================================================
p_bar_group <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("一级 CAZy 微生物贡献（每组内按 ", target_genus, " 比例降序）"),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  ) &
  theme(legend.position = "right")

print(p_bar_group)
ggsave("Bar_CAZy_Genus_byGroup.pdf", p_bar_group, width = 16, height = 6)
ggsave("Bar_CAZy_Genus_byGroup.png", p_bar_group, width = 16, height = 6, dpi = 300)