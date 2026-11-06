library(ComplexUpset)

# 只保留集合列
library(readxl)

data_mat <- read_excel("Cdb.xlsx")
data_mat_sub <- data_mat[, -1]  # 假设第一列是 primary_cluster
p <- upset(
  data_mat_sub,
  intersect = colnames(data_mat_sub),
  min_size = 0,
  width_ratio = 0.1,
  height_ratio = 0.7,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      text = list(size = 3)
    )
  ),
  themes = NULL
)

p