�?
瀹屽叏鍙傝€冭嚜�??
1.鐢熶俊濯?(浼級浠庨浂寮€濮嬪杞綍缁勶紙7\8�??
      2.http://www.jianshu.com/p/324aae3d5ea4
###########################
#浣跨敤DESeq2杩涜宸紓鍩哄洜鍒嗘�?
###########################
#瀵煎叆鏁版嵁锛屾瀯寤? DESeq2 #鎵€闇€�?? DESeqDataSet 瀵硅�?

library(DESeq2)
countData <- raw_count_filt[,2:5]
condition <- factor(c("control","KD","KD","control"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
#杩欎竴姝ュ埌涓嬩竴姝ヤ箣闂村彲浠ヨ繃婊ゆ帀涓€浜沴ow count鏁版嵁锛岃妭鐪佸唴瀛橈紝鎻愰珮杩愯閫熷害

nrow(dds)
#[1] 60492
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
[1] 29493
#浣跨敤DESeq杩涜宸紓琛ㄨ揪鍒嗘瀽锛? DESeq鍖呭惈涓夋
dds <- DESeq(dds)
'''
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
'''
#鐢╮esults鑾峰彇缁撴灉�?? results鐨勫弬鏁伴潪甯哥殑澶氾紝杩欓噷涓嶅ソ鍏蜂綋灞曞紑  浣嗘槸浣犱滑浼氳嚜宸辩湅鐨勫�?

res <- results(dds)
#鍙敤mcols鏌ョ湅姣忎竴椤圭粨鏋滅殑鍏蜂綋鍚箟
mcols(res, use.names = TRUE)
'''
DataFrame with 6 rows and 2 columns
                       type                                     description
                <character>                                     <character>
baseMean       intermediate       mean of normalized counts for all samples
log2FoldChange      results log2 fold change (MLE): condition KD vs control
lfcSE               results         standard error: condition KD vs control
stat                results         Wald statistic: condition KD vs control
pvalue              results      Wald test p-value: condition KD vs control
padj                results                            BH adjusted p-values
'''
#鐢╯ummary鐪嬫弿杩版€х殑缁撴灉
summary(res)
'''
out of 29493 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 3110, 11% 
LFC < 0 (down)   : 2108, 7.1% 
outliers [1]     : 0, 0% 
low counts [2]   : 14867, 50% 
(mean count < 19)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
'''
#鐢讳釜MA鍥撅紝杩樿兘鏍囨敞p鍊兼渶灏忕殑鍩哄�?
#娌℃湁缁忚繃 statistical moderation骞崇紦log2 fold changes鐨勬儏鍐?
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#缁忚繃lfcShrink 鏀剁缉log2 fold change�?? 缁撴灉浼氬ソ鐪嬪緢澶?
res.shrink <- lfcShrink(dds, contrast = c("condition","KD","control"), res=res)
plotMA(res.shrink, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#鍏堟妸宸紓琛ㄨ揪鐨勫熀鍥犳壘鍑烘潵
res.deseq2 <- subset(res, padj < 0.05)
summary(res.deseq2)
'''
out of 4276 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 2570, 60% 
LFC < 0 (down)   : 1706, 40% 
outliers [1]     : 0, 0% 
low counts [2]   : 0, 0% 
(mean count < 19)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
'''
#############################
#浣跨敤edgeR杩涜宸紓鍩哄洜鍒嗘�?
#############################
#绗竴姝ワ細 鏋勫缓DGEList瀵硅�?

library(edgeR)
group <- factor(c("control","KD","KD","control"))
genelist <- DGEList(counts=raw_count_filt[,2:5], group = group)
#绗簩姝ワ細 杩囨�? low counts鏁版嵁銆?
# 绠€鍗曠矖鏆寸殑鏂规�?
keep <- rowSums(genelist$count) > 50
# 鍒╃敤CPM鏍囧噯鍖?
keep <- rowSums(cpm(genelist) > 0.5 ) >=2
table(keep)
'''
keep
FALSE  TRUE 
43944 16548 
'''
genelist.filted <- genelist[keep, ,keep.lib.sizes=FALSE]
summary(genelist.filted)
'''
        Length Class      Mode   
counts  66192  -none-     numeric
samples     3  data.frame list 
'''
#绗笁姝ワ細 鏍规嵁缁勬垚鍋忓�?(composition bias)鏍囧噯鍖栥€俥dgeR鐨刢alcNormFactors鍑芥暟浣跨敤TMM绠楁硶瀵笵GEList鏍囧噯鍖?
genelist.norm <- calcNormFactors(genelist.filted)
#绗洓姝ワ細 瀹為獙璁捐鐭╅樀(Design matrix)�?? 绫讳技浜嶥ESeq2涓殑design鍙傛�?
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
'''
  control KD
1       1  0
2       0  1
3       0  1
4       1  0
attr(,"assign")
[1] 1 1
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"
'''
#绗簲姝ワ細 浼拌绂绘暎鍊硷紙Dispersion锛夈�??
#install.packages(statmod)
genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)
#杩涗竴姝ラ€氳繃quasi-likelihood (QL)鎷熷悎NB妯″瀷锛岀敤浜庤В閲婄敓鐗╁鍜屾妧鏈€у鑷寸殑鍩哄洜鐗瑰紓鎬у彉寮?
fit <- glmQLFit(genelist.Disp, design, robust=TRUE)
head(fit$coefficients)
'''
                  control         KD
ENSG00000000003 -10.194651  -9.981594
ENSG00000000419  -9.537979  -9.436034
ENSG00000000457 -11.244295 -11.521055
ENSG00000000460 -10.456470 -10.345220
ENSG00000001036  -9.600292  -9.267762
ENSG00000001084  -9.706178  -9.399947
'''
#绗叚姝?: 宸紓琛ㄨ揪妫€楠岋�?1锛夈€傝繖涓€姝ヤ富瑕佹瀯寤烘瘮杈冪煩闃?
cntr.vs.KD <- makeContrasts(control-KD, levels=design)
res <- glmQLFTest(fit, contrast=cntr.vs.KD)
ig.edger <- res$table[p.adjust(res$table$PValue, method = "BH") < 0.01, ]
#鍚庣画灏辨槸鎻愬彇鏄捐憲鎬у樊寮傜殑鍩哄洜鐢ㄤ綔涓嬫父鍒嗘瀽锛屽仛涓€浜涘浘鐪嬬湅

topTags(res,n=10)
'''
Coefficient:  1*control -1*KD 
                    logFC   logCPM         F       PValue          FDR
ENSG00000092067  7.318378 5.224942 1918.1434 3.812446e-13 6.308835e-09
ENSG00000175928  5.529875 4.462782 1493.2747 1.392974e-12 1.152547e-08
ENSG00000280322  4.351239 4.750425 1283.8540 3.041784e-12 1.677848e-08
ENSG00000104332  3.458766 4.031653  563.5671 2.102783e-10 8.699212e-07
ENSG00000150893  2.103118 5.219683  501.3174 3.825353e-10 1.266039e-06
ENSG00000036448  1.835316 4.737906  387.6063 1.418224e-09 3.911461e-06
ENSG00000229807  2.726486 9.684564  345.7889 2.530529e-09 5.982171e-06
ENSG00000169136 -1.658243 5.724014  330.3849 3.186983e-09 5.991361e-06
ENSG00000142530 -1.921870 6.145287  328.9378 3.258536e-09 5.991361e-06
ENSG00000062038  2.460219 3.545624  314.3368 4.099064e-09 6.783131e-06
'''
is.de <- decideTestsDGE(res)
summary(is.de)
'''
   1*control -1*KD
-1            3111
0            10169
1             3268
'''
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
#绗叚姝ワ細宸紓琛ㄨ揪妫€楠岋�?2锛夈€傛壘琛ㄨ揪閲忓彉鍖栨瘮杈冨ぇ鐨勫熀鍥狅紝瀵瑰簲鐨勫嚱鏁版�? glmTreat
#B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
tr <- glmTreat(fit, contrast=cntr.vs.KD, lfc=log2(1.5))
#tr <- glmTreat(fit, coef=2, lfc=log2(1.5))
plotMD(tr,status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")

############################
#浣跨敤limma杩涜宸紓鍒嗘�?
############################
#鏁版嵁棰勫鐞嗭�? Limma浣跨敤edgeR鐨凞GEList瀵硅薄锛屽苟涓旇繃婊ゆ柟娉曢兘鏄竴鑷寸殑锛屽搴攅dgeR鐨勭涓€姝?,绗簩姝ワ紝 绗笁姝?

library(edgeR)
library(limma)
group <- factor(c("control","KD","KD","control"))
genelist <- DGEList(counts=raw_count_filt[,2:5], group = group)
### filter base  use CPM
keep <- rowSums(cpm(genelist) > 0.5 ) >=2
table(keep)
'''
keep
FALSE  TRUE 
43944 16548 
'''
genelist.filted <- genelist[keep, ,keep.lib.sizes=FALSE]
### normalizaition
x <- calcNormFactors(genelist.filted, method = "TMM")
#宸紓琛ㄨ揪鍒嗘�?: 浣跨敤鈥漧imma-trend�??

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))
'''
                   logFC  AveExpr        t      P.Value    adj.P.Val        B
ENSG00000112378 7.908618 7.423498 107.1885 3.275833e-12 5.883721e-11 18.31544
ENSG00000070669 7.720042 7.311585 106.4612 3.430632e-12 5.883721e-11 18.28501
ENSG00000143222 7.144495 6.617131 106.3505 3.454928e-12 5.883721e-11 18.28034
ENSG00000105063 8.155105 7.754410 106.0044 3.532139e-12 5.883721e-11 18.26570
ENSG00000167513 8.504830 8.022979 105.9841 3.536740e-12 5.883721e-11 18.26484
ENSG00000204673 7.195007 6.784525 105.7979 3.579168e-12 5.883721e-11 18.25692
ENSG00000161671 7.710280 7.164736 105.7947 3.579908e-12 5.883721e-11 18.25678
ENSG00000088038 7.420823 7.047706 105.6121 3.622078e-12 5.883721e-11 18.24899
ENSG00000126461 7.864620 7.508545 105.4921 3.650112e-12 5.883721e-11 18.24386
ENSG00000057757 7.282897 6.855236 105.3629 3.680589e-12 5.883721e-11 18.23832
'''
#宸紓琛ㄨ揪鍒嗘�?: 浣跨敤鈥漧imma-voom�??

### DGE with voom
v <- voom(genelist.norm, design, plot=TRUE)
#v <- voom(counts, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
all <- topTable(fit, coef=ncol(design), number=10000)
sig.limma <- all[all$adj.P.Val < 0.01, ]
fit <- treat(fit, lfc=log2(1.2))
topTreat(fit, coef=ncol(design))
'''
                   logFC  AveExpr        t      P.Value   adj.P.Val
ENSG00000167658 12.43870 12.37409 151.8180 1.306387e-09 1.40086e-07
ENSG00000080824 12.09251 11.76151 144.6511 1.622197e-09 1.40086e-07
ENSG00000074800 12.40951 11.98962 138.0662 2.007816e-09 1.40086e-07
ENSG00000067225 11.63560 11.39947 136.4730 2.104020e-09 1.40086e-07
ENSG00000184009 11.13266 11.09332 135.8970 2.135837e-09 1.40086e-07
ENSG00000075624 11.36112 11.27684 133.9536 2.284186e-09 1.40086e-07
ENSG00000065978 11.36201 11.03076 133.6619 2.306871e-09 1.40086e-07
ENSG00000109971 11.46223 11.19645 132.2101 2.425811e-09 1.40086e-07
ENSG00000149925 10.80475 10.33477 129.4707 2.652039e-09 1.40086e-07
ENSG00000142676 11.07432 10.81495 128.3201 2.768057e-09 1.40086e-07
'''
#鎻愬彇浜嗗悇绉嶇殑鏄捐憲鎬у熀鍥狅紝姣旇緝灏遍渶瑕佺敤闊︽仼鍥句簡锛屼絾鏄垜鍋忎笉锛屾垜瑕佺敤UpSetR.
#install.packages("UpSetR")
library(UpSetR)
input <- fromList(list(edgeR=rownames(ig.edger), DESeq2=rownames(res.deseq2), limma=rownames(sig.limma)))
summary(input)
'''
edgeR            DESeq2           limma       
Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:1.0000  
Median :0.0000   Median :0.0000   Median :1.0000  
Mean   :0.3173   Mean   :0.3899   Mean   :0.9118  
3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
Max.   :1.0000   Max.   :1.0000   Max.   :1.0000 
'''
upset(input)
#venn�??

vennDiagram(input, include = "both", names = c("edger", "edger", "deseq2"),cex = 1,counts.col = "red")
#limma鍖呬笉鑳界粯鍒跺僵鑹插浘銆傘€傘€傘€傘€傘�??
#install.packages("venneuler")
#install.packages("rJava")
library("rJava")#闇€瑕乯ava......
library(gplots)
library('venneuler')
install.packages("gplots")
venn(input)#涔熶笉鑳芥槸褰╄壊鐨勩€傘€傘€傘€傘€傘€?

library(VennDiagram)
names(input)<-c("edger", "edger", "deseq2")
venn.diagram(input,"VennDiagram.venn.png")
'''
Error in `[[<-.data.frame`(`*tmp*`, i, value = c(1L, 0L)) : 
  replacement has 2 rows, data has 10967
'''




