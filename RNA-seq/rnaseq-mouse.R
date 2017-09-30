#合并表达矩阵

#使用R合并表达矩阵
setwd('G:\\Users\\shenyou\\Desktop\\matrix')

options(stringsAsFactors = FALSE)

#首先将四个文件分别赋值：control1，control2，rep1，rep2

control1 <- read.table("SRR3589959.count", sep = "\t", col.names = c("gene_id", "control1"))

control2 <- read.table("SRR3589961.count", sep= "\t", col.names = c("gene_id", "control2"))

rep1 <- read.table("SRR3589960.count", sep="\t", col.names = c("gene_id", "akap951"))

rep2 <- read.table("SRR3589962.count", sep="\t", col.names = c("gene_id", "akap952"))

#将四个矩阵按照gene_id进行合并，并赋值给raw_count

raw_count <- merge(merge(control1, control2, by="gene_id"), merge(rep1,rep2, by="gene_id"))

#需要将合并的raw_count进行过滤处理，里面有5行需要删除的行，在我们的小鼠的表达矩阵里面，是1,2,48823,48824,48825这5行。并重新赋值给raw_count_filter

raw_count_filt <- raw_count[-48823:-48825, ]

raw_count_filter <- raw_count_filt[-1:-2, ]

#因为我们无法在EBI数据库上直接搜索找到ENSMUSG00000024045.5这样的基因，只能是ENSMUSG00000024045的整数，没有小数点，所以需要进一步替换为整数的形式。

#第一步将匹配到的.以及后面的数字连续匹配并替换为空，并赋值给ENSEMBL

ENSEMBL <- gsub("\\.\\d*", "", raw_count_filter$gene_id)

#将ENSEMBL重新添加到raw_count_filter矩阵

row.names(raw_count_filter) <- ENSEMBL

#看一些基因的表达情况，在UniProt数据库找到AKAP95的id，并从矩阵中找到访问，并赋值给AKAP95变量

AKAP95 <- raw_count_filter[rownames(raw_count_filter)=="ENSMUSG00000024045",]

#查看AKAP95

AKAP95
'''
gene_id control1 control2 akap951 akap952
ENSMUSG00000024045 ENSMUSG00000024045.5     1912     2722    6802    8148
'''
raw_count_filter <- raw_count_filter[ ,-1]
#2.构建dds对象
# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
condition <- factor(c(rep("control",2),rep("akap95",2)), levels = c("control","akap95"))
# 获取count数据
countData <- raw_count_filter[,1:4]
colData <- data.frame(row.names=colnames(raw_count_filter), condition)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
# 查看一下dds的内容
head(dds)
'''
class: DESeqDataSet 
dim: 6 4 
metadata(1): version
assays(1): counts
rownames(6): __no_feature __not_aligned ... ENSMUSG00000000003
ENSMUSG00000000028
rowData names(0):
colnames(4): control1 control2 akap951 akap952
colData names(1): condition
'''
#3.DESeq标准化dds

# normalize 数据
dds2 <- DESeq(dds)
# 查看结果的名称，本次实验中是 "Intercept"，"condition_akap95_vs_control"
resultsNames(dds2)
#[1] "Intercept"                   "condition_akap95_vs_control"
# 将结果用results()函数来获取，赋值给res变量
res <- results(dds2)
# summary一下，看一下结果的概要信息
summary(res)
'''

out of 29252 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 806, 2.8% 
LFC < 0 (down)   : 672, 2.3% 
outliers [1]     : 0, 0% 
low counts [2]   : 13689, 47% 
(mean count < 22)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
'''
#4.4.提取差异分析结果

# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
table(res$padj<0.05)
'''
FALSE  TRUE 
14542  1021
'''
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
# 得到csv格式的差异表达分析结果
write.csv(resdata,file= "control_vs_akap95.cvs",row.names = F)
View(resdata)
library(clusterProfiler)

# 我们是小鼠数据，所以直接安装载入就可以了，当然人类的也是一样。
# 人类的注释数据
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
# 小鼠的注释数据
biocLite("AnnotationDbi")
library(AnnotationDbi)
biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
library('rlang')
install.packages('rlang')
# 这个包应该是clusterProfiler自带的，可以直接载入
library(AnnotationHub)
hub <- AnnotationHub()
keytypes(org.Mm.eg.db)
'''
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
 [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
[9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
[13] "IPI"          "MGI"          "ONTOLOGY"     "ONTOLOGYALL" 
[17] "PATH"         "PFAM"         "PMID"         "PROSITE"     
[21] "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"  
'''
# 进行go分析
ego <- enrichGO(
  gene = row.names(diff_gene_deseq2),
  OrgDb = org.Mm.eg.db,
  keytype = "ENSEMBL",
  ont = "MF"
)
# 气泡图
dotplot(ego, font.size=5)
# 网络图
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
# GO图需要安装额外的包
biocLite("topGO")
biocLite("Rgraphviz")
require(Rgraphviz)
plotGOgraph(ego)

作者：lxmic
链接：http://www.jianshu.com/p/4910d7cec5c8
來源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。