#富集分析
#Y叔的clusterProfiler
#数据筛选，根据padj < 0.05 且Log2FoldChange的绝对值大于1的标准。
library(DESeq2)
#summary(res)
#mcols(res, use.names = TRUE)
deseq2.sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
summary(deseq2.sig)
#安装包下载注释数据（rog.HS.eg.db)
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("clusterProfiler")
biocLite("topGO")
library("topGO")
biocLite("AnnotationHub")
biocLite("Rgraphviz")
library(Rgraphviz)
biocLite("graph")
library(graph)
library(AnnotationHub)
library(clusterProfiler)
ah <- AnnotationHub()
org.hs <- ah[['AH53766']]
###GO富集和GESA
'''
enrichGO(gene=deseq2.sig , OrgDb=org.hs, keytype = "ENTREZID", ont = "MF",
         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, qvalueCutoff = 0.2,
         minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
'''
ego <- enrichGO(
  gene = row.names(deseq2.sig),
  OrgDb = org.hs,
  keytype = "ENSEMBL",
  ont = "MF"
)
#可视化分为冒泡图和网络图等，画起来也就几行代码的事情

dotplot(ego,font.size=5)
enrichMap(ego, vertex.label.cex=0.6, layout=igraph::layout.kamada.kawai)
plotGOgraph(ego)
