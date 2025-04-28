```r
# 加载必要的R包
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(topGO)
library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(topGO)
library(ggplot2)
library(data.table)

# ========== Fig.1a ========== 
# （注：这里重复写了几次 Fig.1a，保留原文）

# 读入差异分析结果
TRMA_VS_TRM <- read.csv("./PRMA_VS_PRM_TrudyOliver.csv")

# 读入GSEA基因集
ALL_GSEA_GMT <- read.gmt("./msigdb.v7.1.symbols.gmt")
colnames(ALL_GSEA_GMT) <- c("term","gene")
ALL_GSEA_GMT$term <- as.character(ALL_GSEA_GMT$term)

# 将人类基因名转换为小鼠基因名
ALL_GSEA_GMT_mus <- ALL_GSEA_GMT
ALL_GSEA_GMT_mus$term <- as.character(ALL_GSEA_GMT_mus$term)
ALL_GSEA_GMT_mus1 <- ALL_GSEA_GMT_mus %>% mutate(mouse_gene = convert_human_to_mouse_symbols(gene))
ALL_GSEA_GMT_mus2 <- ALL_GSEA_GMT_mus1[, c("term", "mouse_gene")]
colnames(ALL_GSEA_GMT_mus2) <- c("term", "gene")
ALL_GSEA_GMT_mus2 <- na.omit(ALL_GSEA_GMT_mus2)

# 准备基因列表
TRMA_VS_TRM1 <- TRMA_VS_TRM[order(TRMA_VS_TRM$log2FoldChange, decreasing=TRUE), ]
TRMA_VS_TRM1$X <- TRMA_VS_TRM1$gene_name
TRMA_VS_TRM_SCLC_vs_Normal_res <- TRMA_VS_TRM1[, c("log2FoldChange", "X")]
TRMA_VS_TRM_SCLC_vs_Normal_res <- na.omit(TRMA_VS_TRM_SCLC_vs_Normal_res)
TRMA_VS_TRM_SCLC_vs_Normal_res_genelist <- TRMA_VS_TRM_SCLC_vs_Normal_res$log2FoldChange
names(TRMA_VS_TRM_SCLC_vs_Normal_res_genelist) <- TRMA_VS_TRM_SCLC_vs_Normal_res$X
TRMA_VS_TRM_SCLC_vs_Normal_res_genelist <- TRMA_VS_TRM_SCLC_vs_Normal_res_genelist[TRMA_VS_TRM_SCLC_vs_Normal_res_genelist != 0]

# GSEA分析
TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA <- GSEA(
  TRMA_VS_TRM_SCLC_vs_Normal_res_genelist, 
  TERM2GENE = ALL_GSEA_GMT_mus2, 
  verbose = TRUE,
  minGSSize = 0,
  pvalueCutoff = 1
)

TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA <- TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA[order(TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA$NES, decreasing = TRUE), ]

# 提取Splicing相关通路
Splicing_pathway <- TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA$ID[grep("SPLICING", TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA$ID)]
TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing <- subset(TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA, ID %in% Splicing_pathway)
TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1 <- TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing[order(TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing$NES, decreasing = FALSE), ][1:6, ]
TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1 <- TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1[order(TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1$NES, decreasing = TRUE), ]
TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1$color <- "Splicing"

# 绘制柱状图
library(ggpubr)
bb <- jdb_palette("brewer_celsius")
color1 <- adjustcolor(bb[8], alpha.f = 0.8)  # 设置颜色透明度60%

p1 <- ggbarplot(
  TRMA_VS_TRM_SCLC_vs_Normal_res_GSEA_splicing1, 
  x = "ID", 
  y = -log10(p.adjust), 
  fill = "color", 
  color = "black",
  title = "TRMA_VS_TRM", 
  outlier.shape = NA, 
  rotate = TRUE,
  lab.size = 1
) +
rotate_x_text(angle = 90) + 
scale_fill_manual(values = color1)

![Fig.1a](./Figures/Fig1/Fig.1a.png)
