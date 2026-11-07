# ============================================================
# 0. 加载包 + 工具函数
# ============================================================
library(dplyr)
library(tidyr)
library(tibble)
library(rstatix)
library(pheatmap)
library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)

strip_path <- function(x) {
  x <- sub("\\s*\\[PATH:.*\\]$", "", x)
  trimws(x)
}

# ============================================================
# 1. 读取 KO 索引，提取脂质代谢 level 3 通路名
# ============================================================
ko_index <- read.delim("KO1-4.txt", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)

lipid_idx_full  <- ko_index %>%
  filter(PathwayL1 == "Metabolism",
         grepl("Lipid metabolism", PathwayL2, ignore.case = TRUE)) %>%
  pull(Pathway) %>% unique()
lipid_idx_short <- strip_path(lipid_idx_full)

cat("索引中脂质代谢 level 3 通路数：", length(lipid_idx_full), "\n")
if (length(lipid_idx_full) == 0) {
  stop("索引里没匹配到 Lipid metabolism，请检查 PathwayL1/PathwayL2 写法")
}

# ============================================================
# 2. 读取 Pathway 丰度表，筛选脂质代谢行
# ============================================================
dat <- read.delim("all_KEGG.Pathway.raw.txt",
                  header = TRUE, check.names = FALSE)
rownames(dat) <- dat[[1]]
dat <- dat[, -1, drop = FALSE]

main_full  <- rownames(dat)
main_short <- strip_path(main_full)

keep <- main_full %in% lipid_idx_full
match_mode <- "full"
if (sum(keep) == 0) {
  keep <- main_short %in% lipid_idx_short
  match_mode <- "short"
}
cat("匹配模式：", match_mode, "｜匹配行数：", sum(keep), "\n")
if (sum(keep) == 0) {
  cat("\n主数据行名示例：\n");     print(head(main_full, 8))
  cat("\n索引脂质通路名示例：\n"); print(head(lipid_idx_full, 8))
  stop("主数据没有任何行匹配脂质代谢通路，请核对命名差异")
}

dat_lipid <- dat[keep, , drop = FALSE]
rownames(dat_lipid) <- main_short[keep]

pathway_lipid <- rownames(dat_lipid)
cat("脂质代谢通路条数：", nrow(dat_lipid), "\n")

# ============================================================
# 3. log2 转换 + 分组信息
# ============================================================
dat_log <- log2(dat_lipid + 1)

group_info <- read.delim("metadata.txt", header = TRUE,
                         stringsAsFactors = FALSE)
desired_order <- c("WM", "NM", "BM")
group_info$Group <- factor(group_info$Group, levels = desired_order)

# ============================================================
# 4. 长表
# ============================================================
long_dat <- dat_log %>%
  as.data.frame() %>%
  mutate(Pathway = rownames(.)) %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "Abundance") %>%
  left_join(group_info, by = c("Sample" = "Samples")) %>%
  filter(!is.na(Group), !is.na(Abundance))

# ============================================================
# 5. 过滤无效通路（至少 2 组有信号）
# ============================================================
long_dat <- long_dat %>%
  group_by(Pathway) %>%
  filter(n_distinct(Group[Abundance > 0]) >= 2) %>%
  ungroup()
cat("过滤后有效脂质通路数：", n_distinct(long_dat$Pathway), "\n")
if (n_distinct(long_dat$Pathway) == 0) {
  stop("所有脂质通路过滤后都不满足条件")
}

# ============================================================
# 6. Kruskal-Wallis + FDR
# ============================================================
kw_res <- long_dat %>%
  group_by(Pathway) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"),
         sig = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ ""
         ))

# ============================================================
# 7. 两两 Wilcoxon + FDR
# ============================================================
pairwise_res <- long_dat %>%
  group_by(Pathway) %>%
  pairwise_wilcox_test(Abundance ~ Group, p.adjust.method = "BH") %>%
  ungroup() %>%
  mutate(sig = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

write.csv(kw_res,       "Lipid_KEGG_KruskalWallis_FDR.csv",   row.names = FALSE)
write.csv(pairwise_res, "Lipid_KEGG_pairwise_Wilcoxon_FDR.csv", row.names = FALSE)

# ============================================================
# 8. 组均值矩阵（列顺序 WM / NM / BM）
# ============================================================
mat_abun <- long_dat %>%
  group_by(Pathway, Group) %>%
  summarise(mean = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean)

present_groups <- intersect(desired_order, colnames(mat_abun))
mat_abun <- mat_abun %>%
  column_to_rownames("Pathway") %>%
  as.matrix()
mat_abun <- mat_abun[, present_groups, drop = FALSE]
mat_abun[is.na(mat_abun)] <- 0

# ============================================================
# 9. 热图：全部脂质通路
# ============================================================
n_show <- nrow(mat_abun)             # 全部；改成 min(30, nrow(...)) 只画前30
mat_abun_top <- mat_abun[seq_len(n_show), , drop = FALSE]

mat_z <- t(scale(t(mat_abun_top)))
mat_z[is.na(mat_z)] <- 0

kw_small <- kw_res %>% filter(Pathway %in% rownames(mat_z))
row_labels <- paste0(
  rownames(mat_z), " ",
  kw_small$sig[match(rownames(mat_z), kw_small$Pathway)]
)

pheatmap(
  mat_z,
  color = colorRampPalette(c("#71B0FD", "white", "#Fd9c99"))(100),
  cluster_rows = TRUE, cluster_cols = FALSE,
  labels_row = row_labels,
  show_colnames = TRUE,
  fontsize_row = 8, fontsize_col = 10,
  main = "Lipid Metabolism KEGG Pathways (level 3)"
)

# ============================================================
# 10. ⭐ ORF → KO（用 KEGG_ko 列，修复原代码的关键点）
# ============================================================
orf_ko <- fread("all_output", header = TRUE, sep = "\t",
                fill = TRUE, na.strings = "")
cat("all_output 列名：\n"); print(names(orf_ko))
# 应有 ORF 列 + KEGG_ko 列

# ============================================================
# 10. ORF → KO（eggNOG-mapper 输出：query + KEGG_ko）
# ============================================================
orf_ko <- fread("all_output", header = TRUE, sep = "\t",
                fill = TRUE, na.strings = "")

# eggNOG-mapper 第一列叫 query，统一改名 ORF
setnames(orf_ko, old = "query", new = "ORF", skip_absent = TRUE)

cat("all_output 列名：\n"); print(names(orf_ko))

orf_ko_long <- orf_ko %>%
  select(ORF, KEGG_ko) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "", KEGG_ko != "-") %>%
  mutate(
    ORF     = trimws(as.character(ORF)),
    KEGG_ko = trimws(as.character(KEGG_ko)),
    KEGG_ko = sub("^ko:", "", KEGG_ko, ignore.case = TRUE)   # 去 ko: 前缀
  ) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = trimws(KEGG_ko)) %>%
  filter(grepl("^K\\d{5}$", KEGG_ko)) %>%
  distinct(ORF, KEGG_ko) %>%
  rename(KO = KEGG_ko)

cat("ORF-KO 映射条数：", nrow(orf_ko_long), "\n")
cat("涉及 ORF 数：",     n_distinct(orf_ko_long$ORF), "\n")
cat("涉及 KO 数：",      n_distinct(orf_ko_long$KO),  "\n")
cat("KO 样例：",  paste(head(unique(orf_ko_long$KO), 5),  collapse = ", "), "\n")
cat("ORF 样例：", paste(head(unique(orf_ko_long$ORF), 5), collapse = ", "), "\n")
# ============================================================
# 11. KO → Pathway（只保留脂质代谢）
# ============================================================
ko_path <- ko_index %>%
  select(KO, Pathway) %>%
  mutate(Pathway = strip_path(Pathway)) %>%
  filter(Pathway %in% pathway_lipid) %>%
  distinct(KO, Pathway)

cat("KO-脂质通路 映射条数：", nrow(ko_path), "\n")
cat("覆盖的脂质通路数：",     n_distinct(ko_path$Pathway), "\n")

# ============================================================
# 12. ORF → Pathway（脂质代谢）
# ============================================================
orf_path_long <- orf_ko_long %>%
  inner_join(ko_path, by = "KO") %>%
  distinct(ORF, Pathway)

cat("ORF-脂质通路 映射条数：", nrow(orf_path_long), "\n")
cat("涉及的 ORF 数：",         n_distinct(orf_path_long$ORF), "\n")

# ============================================================
# 13. ORF 丰度
# ============================================================
ORF_ab <- fread("all_gene.TPM", check.names = FALSE)
ORF_ab_long <- melt(ORF_ab, id.vars = "Name",
                    variable.name = "Samples", value.name = "Abundance")
setnames(ORF_ab_long, "Name", "ORF")
ORF_ab_long <- ORF_ab_long %>% filter(Abundance > 0)

# ============================================================
# 14. ORF 分类
# ============================================================
orf_tax <- fread("all_nucleotide.tax", header = TRUE, sep = "\t",
                 fill = TRUE, na.strings = "")
setnames(orf_tax, "Name", "ORF")

orf_tax_clean <- orf_tax %>%
  filter(!is.na(taxid), taxid != "", taxid != "Unassigned", taxid != "-") %>%
  filter(!Genus %in% c("Homo", "Chordata", "Babesia", "Mus", "Rattus")) %>%
  select(ORF, taxid, Genus, Phylum, Class, Order, Family)

# ============================================================
# 15. 合并（先挂 Pathway 缩行数，再挂分类）
# ============================================================
df <- ORF_ab_long %>%
  inner_join(orf_path_long, by = "ORF") %>%
  inner_join(orf_tax_clean, by = "ORF")

cat("合并后数据行数：", nrow(df), "\n")
if (nrow(df) == 0) {
  cat("TPM∩orf_path_long：", length(intersect(ORF_ab_long$ORF,  orf_path_long$ORF)),  "\n")
  cat("orf_path∩tax     ：", length(intersect(orf_path_long$ORF, orf_tax_clean$ORF)), "\n")
  stop("合并后为空，请检查上面两行 intersect")
}

# ============================================================
# 16. Pathway × Genus 贡献
# ============================================================
path_contrib <- df %>%
  group_by(Pathway, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Pathway) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  arrange(Pathway, desc(Prop)) %>%
  ungroup()

fwrite(as.data.table(path_contrib),
       "Lipid_pathway_genus_contribution_full.tsv",
       sep = "\t", quote = FALSE)

top_path <- path_contrib %>%
  group_by(Pathway) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(Pathway)
cat("脂质通路按总丰度排序：\n"); print(top_path)

# ============================================================
# 17. 颜色方案
# ============================================================
base_colors <- c(
  "Corynebacterium" = "#3188ba",
  "Pseudomonas_E" = "#fbae61",   "Bacteroides" = "#f46d43",
  "Atopostipes" = "#3D9F3C",     "Mammaliicoccus" = "#9d0e43",
  "Burkholderia" = "#d53e4f",    "Dietzia" = "#abd7a2",
  "Enterococcus" = "#A1CDEC",    "Alloiococcus" = "#f1c8a8",
  "Nitriliruptor" = "#e9792e",   "Streptococcus" = "#f5f8be",
  "Shewanella" = "#c7b3d6",      "Staphylococcus" = "#fddbc7",
  "Neisseria" = "#c5d9e8",       "Nostoc" = "#b8e186",
  "Arthrobacter" = "#f4a582",    "Enterobacter" = "#92c5de",
  "Providencia" = "#d6604d",     "Acinetobacter" = "#4393c3","Trueperella"     = "#66c2a5",
  "Pseudomonas"     = "#abdda4",
  "Others" = "#bdbdbd",          "Unassigned" = "#636363"
)

# 鲁棒颜色生成：从 levels 直接构造命名向量，杜绝 name 不匹配
make_colors <- function(lv, palette = base_colors) {
  if (is.null(lv) || length(lv) == 0) return(setNames(character(0), character(0)))
  out  <- setNames(rep(NA_character_, length(lv)), lv)
  hit  <- intersect(lv, names(palette))
  out[hit] <- palette[hit]
  miss <- setdiff(lv, hit)
  if (length(miss) > 0) {
    out[miss] <- colorRampPalette(brewer.pal(12, "Set3"))(length(miss))
  }
  stopifnot(!any(is.na(out)))
  out
}

# ============================================================
# 18. 饼图 1：含 Others + Unassigned
# ============================================================
df_genus_total <- path_contrib %>%
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

top_n_genera <- 10
top_genera <- df_genus_total %>%
  filter(Genus != "Unassigned") %>%
  slice_max(Prop, n = top_n_genera) %>%
  pull(Genus)

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
df_pie_all$Group <- factor(as.character(df_pie_all$Group), levels = group_order)
df_pie_all <- df_pie_all %>%
  mutate(label = paste0(Group, "\n", percent(Prop, accuracy = 0.1)))

pie_colors_all <- make_colors(levels(df_pie_all$Group))

p_pie_all <- ggplot(df_pie_all, aes(x = "", y = Prop, fill = Group)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pie_colors_all) +
  theme_void() +
  labs(fill = "Genus",
       title = "脂质代谢通路 属贡献总览（含 Others / Unassigned）")
print(p_pie_all)
ggsave("Pie_Lipid_Genus_All.pdf", p_pie_all, width = 8, height = 8)
ggsave("Pie_Lipid_Genus_All.png", p_pie_all, width = 8, height = 8, dpi = 300)

# ============================================================
# 19. 饼图 2：去掉 Others + Unassigned
# ============================================================
df_pie_clean <- df_genus_total %>%
  filter(!Genus %in% c("Unassigned", "Others")) %>%
  slice_max(Prop, n = top_n_genera) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  arrange(desc(Prop))
df_pie_clean$Genus <- factor(as.character(df_pie_clean$Genus),
                             levels = df_pie_clean$Genus)
df_pie_clean <- df_pie_clean %>%
  mutate(label = paste0(Genus, "\n", percent(Prop, accuracy = 0.1)))

pie_colors_clean <- make_colors(levels(df_pie_clean$Genus))

p_pie_clean <- ggplot(df_pie_clean, aes(x = "", y = Prop, fill = Genus)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pie_colors_clean) +
  theme_void() +
  labs(fill = "Genus",
       title = "脂质代谢通路 属贡献（去 Others / Unassigned）")
print(p_pie_clean)
ggsave("Pie_Lipid_Genus_Clean.pdf", p_pie_clean, width = 8, height = 8)
ggsave("Pie_Lipid_Genus_Clean.png", p_pie_clean, width = 8, height = 8, dpi = 300)

# ============================================================
# 20. 堆叠柱状图：x = 脂质通路，fill = Genus
# ============================================================
df_clean <- path_contrib %>%
  filter(!Genus %in% c("Others", "Unassigned", "Unclassified"),
         !is.na(Genus), Genus != "") %>%
  filter(Pathway %in% top_path)

top10_genera <- df_clean %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

df_plot <- df_clean %>% filter(Genus %in% top10_genera)

df_plot <- df_plot %>%
  group_by(Pathway) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  ungroup()

# 通路排序（按总丰度降序）
path_order <- df_plot %>%
  group_by(Pathway) %>%
  summarise(t = sum(Abundance), .groups = "drop") %>%
  arrange(desc(t)) %>% pull(Pathway)
df_plot$Pathway <- factor(as.character(df_plot$Pathway), levels = path_order)

# 属排序（按总量降序）
genus_order <- df_plot %>%
  group_by(Genus) %>%
  summarise(total = sum(Prop), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(Genus)
df_plot$Genus <- factor(as.character(df_plot$Genus), levels = genus_order)

mycol <- make_colors(levels(df_plot$Genus))

cat("\n===== 画图前检查 =====\n")
cat("nrow(df_plot)   :", nrow(df_plot),   "\n")
cat("levels(Genus)   :", paste(levels(df_plot$Genus),   collapse = " | "), "\n")
cat("names(mycol)    :", paste(names(mycol),            collapse = " | "), "\n")
cat("共同元素个数    :", length(intersect(names(mycol),
                                    as.character(df_plot$Genus))), "\n")

##########
########## 分组堆叠柱状图（像 CAZy 一样按组分面并独立排序） ##########

library(patchwork) # 确保加载了拼接包

# ============================================================
# 1. 挂载分组信息并计算分组贡献
# ============================================================
# 将 metadata 的分组信息 join 到之前的 df 上
df_grp <- df %>%
  left_join(group_info %>% select(Samples, Group), by = "Samples") %>%
  filter(!is.na(Group))

cat("挂载分组后数据行数：", nrow(df_grp), "\n")

# 按 组、通路、属 汇总丰度并计算组内相对比例
path_contrib_grp <- df_grp %>%
  group_by(Group, Pathway, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Pathway) %>%
  filter(sum(Abundance) > 0) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  ungroup()

# 导出分组数据表
fwrite(as.data.table(path_contrib_grp),
       "Lipid_pathway_genus_contribution_byGroup.tsv",
       sep = "\t", quote = FALSE)

# ============================================================
# 2. 取 Top 10 属（跨所有组合并，保证三组颜色/图例一致）
# ============================================================
df_clean_grp <- path_contrib_grp %>%
  filter(!Genus %in% c("Others", "Unassigned", "Unclassified"),
         !is.na(Genus), Genus != "")

top10_genera_grp <- df_clean_grp %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

df_plot_grp <- df_clean_grp %>% filter(Genus %in% top10_genera_grp)

# 每个组内、每个 Pathway 重新归一化
df_plot_grp <- df_plot_grp %>%
  group_by(Group, Pathway) %>%
  mutate(Prop = Abundance / sum(Abundance)) %>%
  ungroup()

# ============================================================
# 3. 颜色方案
# ============================================================
# 使用你在第 17 步定义的 base_colors 和 make_colors 函数
mycol_grp <- make_colors(unique(as.character(df_plot_grp$Genus)))

# 设定排序目标（棒状杆菌属）
target_genus <- "Corynebacterium"

# ============================================================
# 4. 分组画堆叠图（每组内独立按 target_genus 的 Prop 降序）
# ============================================================
plot_list <- list()

for (g in levels(df_plot_grp$Group)) {
  
  d <- df_plot_grp %>% filter(Group == g)
  if (nrow(d) == 0) next
  
  # 计算当前组内，每条 Pathway 中 target_genus 的相对比例 (Prop)
  if (target_genus %in% unique(as.character(d$Genus))) {
    ord_df <- d %>%
      mutate(Genus_chr = as.character(Genus)) %>%
      filter(Genus_chr == target_genus) %>%
      group_by(Pathway) %>%
      summarise(prop_val = sum(Prop, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(prop_val))
    
    path_order <- c(ord_df$Pathway, setdiff(unique(as.character(d$Pathway)), ord_df$Pathway))
  } else {
    # 兜底：如果没有 Corynebacterium，按总丰度排
    path_order <- d %>%
      group_by(Pathway) %>%
      summarise(t = sum(Abundance), .groups = "drop") %>%
      arrange(desc(t)) %>% pull(Pathway)
  }
  
  # 设定当前组的因子水平
  d$Pathway <- factor(as.character(d$Pathway), levels = path_order)
  
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
  p <- ggplot(d, aes(x = Pathway, y = Prop, fill = Genus)) +
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
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  plot_list[[g]] <- p
}

# ============================================================
# 5. 拼图并保存
# ============================================================
p_bar_group <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("脂质代谢通路 微生物贡献（每组内按 ", target_genus, " 比例降序）"),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  ) &
  theme(legend.position = "right")

print(p_bar_group)
ggsave("Bar_Lipid_Pathway_Genus_byGroup.pdf", p_bar_group, width = 18, height = 6)
ggsave("Bar_Lipid_Pathway_Genus_byGroup.png", p_bar_group, width = 18, height = 6, dpi = 300)

