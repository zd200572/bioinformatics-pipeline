#�ϲ��������

#ʹ��R�ϲ��������
setwd('G:\\Users\\shenyou\\Desktop\\matrix')

options(stringsAsFactors = FALSE)

#���Ƚ��ĸ��ļ��ֱ�ֵ��control1��control2��rep1��rep2

control1 <- read.table("SRR3589959.count", sep = "\t", col.names = c("gene_id", "control1"))

control2 <- read.table("SRR3589961.count", sep= "\t", col.names = c("gene_id", "control2"))

rep1 <- read.table("SRR3589960.count", sep="\t", col.names = c("gene_id", "akap951"))

rep2 <- read.table("SRR3589962.count", sep="\t", col.names = c("gene_id", "akap952"))

#���ĸ�������gene_id���кϲ�������ֵ��raw_count

raw_count <- merge(merge(control1, control2, by="gene_id"), merge(rep1,rep2, by="gene_id"))

#��Ҫ���ϲ���raw_count���й��˴�����������5����Ҫɾ�����У������ǵ�С��ı���������棬��1,2,48823,48824,48825��5�С������¸�ֵ��raw_count_filter

raw_count_filt <- raw_count[-48823:-48825, ]

raw_count_filter <- raw_count_filt[-1:-2, ]

#��Ϊ�����޷���EBI���ݿ���ֱ�������ҵ�ENSMUSG00000024045.5�����Ļ���ֻ����ENSMUSG00000024045��������û��С���㣬������Ҫ��һ���滻Ϊ��������ʽ��

#��һ����ƥ�䵽��.�Լ��������������ƥ�䲢�滻Ϊ�գ�����ֵ��ENSEMBL

ENSEMBL <- gsub("\\.\\d*", "", raw_count_filter$gene_id)

#��ENSEMBL�������ӵ�raw_count_filter����

row.names(raw_count_filter) <- ENSEMBL

#��һЩ����ı����������UniProt���ݿ��ҵ�AKAP95��id�����Ӿ������ҵ����ʣ�����ֵ��AKAP95����

AKAP95 <- raw_count_filter[rownames(raw_count_filter)=="ENSMUSG00000024045",]

#�鿴AKAP95

AKAP95
'''
gene_id control1 control2 akap951 akap952
ENSMUSG00000024045 ENSMUSG00000024045.5     1912     2722    6802    8148
'''
raw_count_filter <- raw_count_filter[ ,-1]
#2.����dds����
# ��һ���ܹؼ���Ҫ����condition���������ӣ������������ƣ�С�������ж�����ʹ����飬�������ظ�
condition <- factor(c(rep("control",2),rep("akap95",2)), levels = c("control","akap95"))
# ��ȡcount����
countData <- raw_count_filter[,1:4]
colData <- data.frame(row.names=colnames(raw_count_filter), condition)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
# �鿴һ��dds������
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
#3.DESeq��׼��dds

# normalize ����
dds2 <- DESeq(dds)
# �鿴��������ƣ�����ʵ������ "Intercept"��"condition_akap95_vs_control"
resultsNames(dds2)
#[1] "Intercept"                   "condition_akap95_vs_control"
# �������results()��������ȡ����ֵ��res����
res <- results(dds2)
# summaryһ�£���һ�½���ĸ�Ҫ��Ϣ
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
#4.4.��ȡ����������

# ��ȡpadj��pֵ��������У��У�����ֵ��С��0.05�����ﱶ��ȡ��2Ϊ���������1����С��-1�Ĳ���������
table(res$padj<0.05)
'''
FALSE  TRUE 
14542  1021
'''
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
# �õ�csv��ʽ�Ĳ������������
write.csv(resdata,file= "control_vs_akap95.cvs",row.names = F)
View(resdata)
library(clusterProfiler)

# ������С�����ݣ�����ֱ�Ӱ�װ����Ϳ����ˣ���Ȼ�����Ҳ��һ����
# �����ע������
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
# С���ע������
biocLite("AnnotationDbi")
library(AnnotationDbi)
biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
library('rlang')
install.packages('rlang')
# �����Ӧ����clusterProfiler�Դ��ģ�����ֱ������
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
# ����go����
ego <- enrichGO(
  gene = row.names(diff_gene_deseq2),
  OrgDb = org.Mm.eg.db,
  keytype = "ENSEMBL",
  ont = "MF"
)
# ����ͼ
dotplot(ego, font.size=5)
# ����ͼ
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
# GOͼ��Ҫ��װ����İ�
biocLite("topGO")
biocLite("Rgraphviz")
require(Rgraphviz)
plotGOgraph(ego)

���ߣ�lxmic
���ӣ�http://www.jianshu.com/p/4910d7cec5c8
��Դ������
����Ȩ���������С���ҵת������ϵ���߻����Ȩ������ҵת����ע��������