# -------------------------
# 1.R包载入
# -------------------------
library(ggplot2)
library(gghalves)      
library(dplyr)
library(ggsignif)      
library(ggpubr)       
library(readxl)
# 2. 读取输入文件
# 读取OTU表
alpha_div <- read.table("bracken.S.alpha", header = TRUE, row.names = 1, sep = "\t", comment.char = "")

# 读取实验设计
metadata = read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

# 3.设置 Group 的顺序（保证 good_musk 在最左边）
metadata$Group <- factor(metadata$Group, levels = c("WM", "NM", "BM"))

#4.
alpha_div <- alpha_div %>% 
  mutate(Sample = rownames(alpha_div))

metadata <- metadata %>% 
  mutate(Sample = rownames(metadata))

# 左连接，保证 alpha_div 的顺序
df <- left_join(alpha_div, metadata, by = "Sample")

# 检查列名
colnames(df)

# -------------------------
# 5. 绘制图形
# -------------------------
p<-ggplot(df,aes(x= Group, y=shannon ,fill=Group,color=Group))+
  scale_fill_manual(values = c("#dfc27d","#bf812d","#543005"))+   
  scale_colour_manual(values = c("#dfc27d","#bf812d","#543005"))
p

p1<-p+geom_half_violin(position=position_nudge(x=0.1,y=0),                       
                       side='R',adjust=1.2,trim=F,color=NA,alpha=0.5)
p1

#添加散点图#调整散点
p2<-p1+geom_half_point(position=position_nudge(x=-0.35,y=0),size =3,                        
                       shape =19,range_scale = 0.5,alpha=0.5)
p2
#添加箱线图
p3<-p2+geom_boxplot(outlier.shape = NA, #隐藏离群点；                    
                    width =0.1,                    
                    alpha=0.5)
p3


windowsFonts(A=windowsFont("Times New Roman"), 
             B=windowsFont("Arial"))
p4 <- p3+theme_bw()+theme(panel.grid=element_blank())+  
  theme(axis.ticks.length=unit(-0.25, "cm"),         
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),         
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )+  
  theme(text=element_text(family="A",size=20))
p4


df$shannon <- as.numeric(df$shannon)

# 添加所有两两组比较的显著性
# -------------------------
comparisons <- list(c("NM", "BM"),
                    c("NM", "WM"),
                    c("BM", "WM"))

# ===============================
# 6. 执行 Wilcoxon 检验并输出结果表
# ===============================
#df$shannon <- as.numeric(df$shannon)
test_results <- pairwise.wilcox.test(df$shannon, df$Group, p.adjust.method = "BH")  # Benjamini-Hochberg 校正

# 将结果整理为数据框便于导出
wilcox_df <- as.data.frame(as.table(test_results$p.value))
wilcox_df <- na.omit(wilcox_df)  # 去掉 NA
wilcox_df <- wilcox_df %>%
  rename(Group1 = Var1, Group2 = Var2, p_value = Freq) %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# 打印结果表
print(wilcox_df)
#write.table(wilcox_df, "alpha_shannon_wilcox_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 调整标记位置
max_val <- max(df$shannon)
#max_val <- max(df$Shannon)
p5 <- p4 + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",
  label = "p.signif",
  label.y = c(max_val * 1.05, max_val * 1.1, max_val * 1.15),
  vjust = 0.8   # 👈 星号往下移（数值越大越往下）
)
# -------------------------
# 绘制
# -------------------------
print(p5)