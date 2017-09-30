# import data if sample are small
options(stringsAsFactors = FALSE)
control <- read.table("H:\\biodata\\SRR3589956.count",
                       sep="\t", col.names = c("gene_id","control"))
rep1 <- read.table("H:\\biodata\\SRR3589957.count",
                    sep="\t", col.names = c("gene_id","rep1"))
rep2 <- read.table("H:\\biodata\\SRR3589958.count",
                    sep="\t",col.names = c("gene_id","rep2"))
# merge data and delete the unuseful row

raw_count <- merge(merge(control, rep1, by="gene_id"), rep2, by="gene_id")
raw_count_filt <- raw_count[-1:-5,]

ENSEMBL <- gsub("(.*?)\\.\\d*?_\\d", "\\1", raw_count_filt$gene_id)
row.names(raw_count_filt) <- ENSEMBL
## the sample problem

delta_mean <- abs(mean(raw_count_filt$rep1) - mean(raw_count_filt$rep2))

sampleNum <- length(raw_count_filt$control)
sampleMean <- mean(raw_count_filt$control)
control2 <- integer(sampleNum)
for (i in 1:sampleNum){  if(raw_count_filt$control[i] < sampleMean){
    control2[i] <- raw_count_filt$control[i] + abs(raw_count_filt$rep1[i] - raw_count_filt$rep2[i])
  }  else{
    control2[i] <- raw_count_filt$control[i] + rpois(1,delta_mean)
  }
}
# add data to raw_count

raw_count_filt$control2 <- control2