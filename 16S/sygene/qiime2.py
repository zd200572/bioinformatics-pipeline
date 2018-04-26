/home/zjd/qiime2/171213_16s
gzip *
source activate qiime2-2018.2
##1## import data
qiime tools import  \
 --type 'SampleData[PairedEndSequencesWithQuality]'  \
 --input-path 171213_16s-manifest  \
 --output-path 171213_16s.qza  \
 --source-format PairedEndFastqManifestPhred33
##2## quality control
#visualization
qiime demux summarize \
--i-data 171213_16s.qza\
 --o-visualization 171213_16s.qzv
#filter
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 171213_16s.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 250 \
--o-representative-sequences rep-seqs-dada2-0f250r.qza \
--o-table table-dada2-0f250r.qza \
--p-n-threads 40
#qzv
qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv
#export
qiime tools export \
  table-dada2.qza \
  --output-dir exported-feature-table
##3##代表序列统计
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv
############建树：用于多样性分析##########
# 多序列比对
qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza
# 移除高变区
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
# 建树
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
# 无根树转换为有根树
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
######Alpha#########
qiime diversity core-metrics \
  --i-table table-dada2.qza \
  --p-sampling-depth 5000 \
  --output-dir core-metrics-results-5000 \
  --m-metadata-file sample-metadata.tsv
# --i-phylogeny rooted-tree.qza \

# 统计evenness组间差异是否显著#######报错metadata缺少分类信息
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
# 下载物种注释
#wget -O "gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2017.7/common/gg-13-8-99-515-806-nb-classifier.qza"
# 物种分类
qiime feature-classifier classify-sklearn \
  --i-classifier ../gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs-dada2-0f250r.qza \
  --o-classification taxonomy-0f250r.qza
# 物种结果转换表格，可用于查看
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
# 物种分类柱状图
qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy-0f250r.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv  
  
source deactivate