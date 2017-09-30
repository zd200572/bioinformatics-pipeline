setwd('G:\\Users\\shenyou\\Desktop\\DEG')
a=read.table('siAkap95.mm10.txt')
library(reshape2)
colnames(a)=c('sample','gene','reads')
exprSet=dcast(a,gene~sample)
exprSet=exprSet[-c(1:6),]
colnames(exprSet)
colnames(exprSet) <- c('gene','control_rep2','siAkap95_rep2','control_rep1','siAkap95_rep1')
write.csv(exprSet,'exprSet.counts.csv',row.names = F)

rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
tmp_exprSet <- exprSet
group_list =c('control','siAkap95','control','siAkap95')
keep <- rowSums(edgeR::cpm(exprSet)>10) >= 1
exprSet <- exprSet[keep,]  ## just keep 10507 genes !!!
#source("http://bioconductor.org/biocLite.R")
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#biocLite("edgeR")
library(edgeR)
d <- DGEList(counts=exprSet,group=factor(group_list))
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))

d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)
write.csv(d1$pseudo.counts,file = 'edgeR.pseudo.counts.csv')

et12 <- exactTest(d1) 

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)


nrDEG=topTags(et12, n=nrow(exprSet))
nrDEG=as.data.frame(nrDEG)
write.csv(nrDEG,"edger_classic.results.csv",quote = F)

suppressMessages(library(DESeq2))
exprSet=ceiling(exprSet)
(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list))
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds) 
plotDispEsts(dds, main="Dispersion plot")

rld <- rlogTransformation(dds)
write.csv(assay(rld) ,file = 'DESeq2.pseudo.counts.csv')

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
write.csv(resOrdered,"deseq2.results.csv",quote = F)

