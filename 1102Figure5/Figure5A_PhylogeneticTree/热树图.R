# ============================================================
# 使用 metacoder 绘制热树图（直接显示在 R 中，带节点标签）
# 数据：tax_self.txt 和 abundance.tsv
# ============================================================
# ============================================================

library(metacoder)

# 1. 读取数据
tax <- read.table("tax_self.txt", header = TRUE, sep = "\t", 
                  stringsAsFactors = FALSE, check.names = FALSE)
abund <- read.table("abundance.tsv", header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

# 2. 拼接 lineage（按 ; 分隔，便于 metacoder 识别层级）
tax$lineage <- apply(tax[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus")], 1,
                     function(x) paste(x, collapse = ";"))

# 3. 构建输入数据框
rownames(tax) <- tax$ID
dat2 <- cbind(otu_id = tax$ID, lineage = tax$lineage, 
              abund[ tax$ID, , drop = FALSE])

# 4. 解析数据
obj <- parse_tax_data(dat2,
                      class_cols = "lineage",
                      class_sep = ";")

# 5. 过滤低丰度（去除因测序错误引入的极低 read 数据）
obj$data$tax_data <- zero_low_counts(obj, dataset = "tax_data", min_count = 5)
no_reads <- rowSums(obj$data$tax_data[, -1]) == 0
obj <- filter_obs(obj, target = "tax_data", !no_reads, drop_taxa = TRUE)

# 6. 抽平并计算分类水平丰度
obj$data$tax_data <- calc_obs_props(obj, "tax_data")

# 获取所有样本的列名
sample_cols <- colnames(abund)

# 计算分类层次的总相对丰度
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data", cols = sample_cols)

# ==================== 核心计算：总体流行率百分比 ====================

# 7. 统计该分类在所有样本中出现的绝对次数（检出样本数）
# 注意：这里去掉 groups 参数，计算所有样本的总体检出率
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", cols = sample_cols)

# 8. 计算总体流行率百分比（检出样本数 / 总样本数 * 100）
total_samples <- ncol(abund) # 总样本数
# 直接在对象中添加一列 prevalence_pct，数值范围 0~100
obj$data$tax_occ$prevalence_pct <- (obj$data$tax_occ$n_samples / total_samples) * 100

# ==================== 定制节点标签（去除前缀） ====================
raw_names <- obj$taxon_names()
clean_names <- gsub("^[a-z]__", "", raw_names)
clean_names[clean_names == ""] <- NA
obj$set_taxon_names(clean_names)

# ==================== 绘图：节点颜色=流行率，节点大小=OTU数 ====================
set.seed(1)
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,                              # 节点大小 = OTU 物种数量
          node_size_trans = "log10 area",                 # 🔴 核心修改：使用 log10 面积转换
          node_color = obj$data$tax_occ$prevalence_pct,   # 节点颜色 = 总体流行率 (%)
          
          node_size_axis_label = "Species count (log10)", # 图例说明加上 log10
          node_color_axis_label = "Samples with reads (%)", 
          
          # 圆形布局参数
          layout = "reingold-tilford", 
          aspect_ratio = 1,                        
          
          # 美观防重叠参数
          node_size_range = c(0.005, 0.06),
          node_label_size_range = c(0.03, 0.03), 
          overlap_avoidance = 0.8,
          repel_labels = TRUE,
          repel_force = 3,
          repel_iter = 8000,
          node_label_max = 200)
heat_tree(obj, node_label = taxon_names, node_size = n_obs, node_color = obj$data$tax_occ$prevalence_pct,
          node_size_trans = "log10 area")
