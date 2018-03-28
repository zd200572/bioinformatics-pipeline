#! /bin/bash
#############################################################
# Title: Pipeline of 16S amplicon by Hiseq2500-PE250
# Author: Yong-Xin Liu
# Wechat: yongxinliu
# E-mail: yxliu@genetics.ac.cn
# Website: http://bailab.genetics.ac.cn/
# Date: 3/16/2017
# Version: 1.3
# Enviroment: Ubuntu 16.04x64, qiime 1.9.1, usearch8
# Description: Script for automatic from clean data and mapping file to otu table, taxonomy, alpha and beta raw result
# Requre file list:
# 1. clean data: clean_data/PE250_1/2.fq.gz
# 2. mapping file: doc/PE250.mappingfile.txt
#############################################################
# Hiseq2500 PE250 of bacterial 16S, data in clean_data/ and each library mapping file in doc/ 

# 建立小数据集用于测序分析
cd ~/test/example_PE250
cp ~/ath/jt.terpene.16S/clean_data/T131_*.gz ./
less T131_1.fq.gz | wc -l # 32637276 总共有32637276条序列，了10000000条作测试
zless T131_1.fq.gz|head -n 10000000 |gzip > PE250_1.fq.gz 
zless T131_2.fq.gz|head -n 10000000 |gzip > PE250_2.fq.gz


# 1. 下载实验数据双端测序数据PE250_1/2.fq.gz，和实验设计mappingfile.txt至当前工作目录
mkdir temp result

# 2. 验证实验设计是否有错误
validate_mapping_file.py -m mappingfile.txt

# 3. 合并双端数据
join_paired_ends.py -f PE250_1.fq.gz -r PE250_2.fq.gz -m fastq-join -o temp/PE250_join

# 4. 提取barcode
extract_barcodes.py -f temp/PE250_join/fastqjoin.join.fastq \
	-m mappingfile.txt \
	-o temp/PE250_barcode \
	-c barcode_paired_stitched --bc1_len 0 --bc2_len 6 -a --rev_comp_bc2
# 5. 拆分样品
split_libraries_fastq.py -i temp/PE250_barcode/reads.fastq \
	-b temp/PE250_barcode/barcodes.fastq \
	-m mappingfile.txt \
	-o temp/PE250_split/ \
	-q 20 --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type 6

# 6. Remove adaptor
cutadapt -g AACMGGATTAGATACCCKG -a GGAAGGTGGGGATGACGT -e 0.15 -m 300 --discard-untrimmed temp/PE250_split/seqs.fna -o temp/PE250_P5.fa

# 7. format to Usearch
sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' temp/PE250_P5.fa > temp/seqs_usearch.fa

# 8. Dereplication
./usearch10 -derep_fulllength temp/seqs_usearch.fa \
	-fastaout temp/seqs_unique.fa \
	-minuniquesize 2 -sizeout

# 9. Cluster OTU
./usearch10 -cluster_otus temp/seqs_unique.fa -otus temp/otus.fa -uparseout temp/uparse.txt -relabel Otu

# 10. Remove chimeras
./usearch10 -uchime2_ref temp/otus.fa \
	-db rdp_gold.fa \
	-chimeras temp/otus_chimeras.fa \
	-notmatched temp/otus_rdp.fa \
	-uchimeout temp/otus_rdp.uchime \
	-strand plus -mode sensitive -threads 96

grep '>' temp/otus_chimeras.fa|sed 's/>//g' > temp/otus_chimeras.id

filter_fasta.py -f temp/otus.fa -o temp/otus_non_chimera.fa -s temp/otus_chimeras.id -n
grep '>' -c temp/otus_non_chimera.fa # 2835

# 11. 去除非细菌序列
# http://greengenes.secondgenome.com/downloads/database/13_5
# gunzip gg_13_5_pynast.fasta.gz # 原文件500M，解压后有9个G，这么大的数据，怎么比较
#nohup align_seqs.py -i temp/otus_non_chimera.fa -t gg_13_5_pynast.fasta -o temp/aligned/ &bg # 不超过3h完成
#grep '>' temp/aligned/otus_non_chimera_failures.fasta # 1844 又三分之二不是细菌吗？

# QIIME提供的下载链接
# 下载Greengene最新数据库，320MB
wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# 解压数据包后大小3.4G
tar xvzf gg_13_8_otus.tar.gz
# 将OTU与97%相似聚类的代表性序列多序列比对，大约8min
time align_seqs.py -i temp/otus_non_chimera.fa -t gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
# 无法比对细菌的数量
grep -c '>' temp/aligned/otus_non_chimera_failures.fasta # 1860
# 获得不像细菌的OTU ID
grep '>' temp/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_non_chimera_failures.id
# 过滤非细菌序列
filter_fasta.py -f temp/otus_non_chimera.fa -o temp/otus_rdp_align.fa -s temp/aligned/otus_non_chimera_failures.id -n
# 看我们现在还有多少OTU:975
grep '>' -c temp/otus_rdp_align.fa 


# 12. Generate representitive sequences(rep seqs) and OTU table, remove low abundance samples
# 重命名OTU，这就是最终版的代表性序列，即Reference
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
# 生成OTU表
./usearch10 -usearch_global temp/seqs_usearch.fa -db result/rep_seqs.fa -otutabout temp/otu_table.txt -biomout temp/otu_table.biom -strand plus -id 0.97 -threads 10
# 结果汇总 01:20 141Mb   100.0% Searching seqs_usearch.fa, 32.3% matched
# 默认10线程，用时1分20秒，有32.3%的序列匹配到OTU上；用30线程反而用时3分04秒


# 13. 物种注释Taxonomy assignment
assign_taxonomy.py -i result/rep_seqs.fa \
	-r gg_13_8_otus/rep_set/97_otus.fasta \
	-t gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
	-m rdp -o result

# 14. OTU格式转换，添加信息

biom convert -i temp/otu_table.txt \
	-o result/otu_table.biom \
	--table-type="OTU table" --to-json

biom add-metadata -i result/otu_table.biom \
	--observation-metadata-fp result/rep_seqs_tax_assignments.txt \
	-o result/otu_table_tax.biom \
	--sc-separated taxonomy --observation-header OTUID,taxonomy 

# 15. OTU表数据筛选

# 按样品数据量过滤：选择counts>3000的样品
filter_samples_from_otu_table.py -i result/otu_table_tax.biom -o result/otu_table2.biom -n 3000
# 查看过滤后结果：只有25个样品，975个OTU
biom summarize-table -i result/otu_table2.biom

# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i result/otu_table2.biom -o result/otu_table3.biom
# 查看过滤后结果：只有25个样品，346个OTU
biom summarize-table -i result/otu_table3.biom

# 按物种筛选OTU表：去除p__Chloroflexi菌门
filter_taxa_from_otu_table.py -i result/otu_table3.biom -o result/otu_table4.biom -n p__Chloroflexi
# 查看过滤后结果：只有25个样品，307个OTU
biom summarize-table -i result/otu_table4.biom

# 转换最终biom格式OTU表为文本OTU表格
biom convert -i result/otu_table4.biom -o result/otu_table4.txt --table-type="OTU table" --to-tsv
# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt
# 筛选最终OTU表中对应的OTU序列
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep4.fa -b result/otu_table4.biom


# 16. 进化
clustalo -i result/rep_seqs4.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=30
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree # generate tree by FastTree

# 17. Alpha多样性
biom summarize-table -i result/otu_table4.biom
single_rarefaction.py -i result/otu_table4.biom -o temp/otu_table_rare.biom -d 2797
alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 18. Beta diversity
normalize_table.py -i result/otu_table4.biom -o temp/otu_table_css.biom -a CSS
biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
sed -i 's/^\t//g' result/beta/*


# 19. 按物种分类级别分类汇总


summarize_taxa.py -i result/otu_table4.biom -o result/sum_taxa # summary each level percentage
rm result/sum_taxa/*.biom
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* # format for R read



# 20. 筛选可展示的进化树

filter_otus_from_otu_table.py --min_count_fraction 0.001 -i result/otu_table4.biom -o temp/otu_table_k1.biom
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/otu_table_k1.biom 
grep -c '>' temp/tax_rep_seqs.fa # 104
clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=30
make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/; /\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax


# 21. 其它

# 将mappingfile转换为R可读的实验设计
sed 's/#//' mappingfile.txt > result/design.txt
# 转换文本otu_table格式为R可读
sed '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt > result/otu_table.txt
# 转换物种注释信息为制表符分隔，方便R读取
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt


# 0. Set enviroment and format sample name
wd=/mnt/bai/yongxin/ath/jt.terpene.16S/v1.3 # working directory
cd $wd
# For practice, quickly prepare data and mapping file in your own directory
#ln -s /mnt/bai/yongxin/ath/jt.terpene.16S/clean_data/
#ln -s /mnt/bai/yongxin/ath/jt.terpene.16S/doc/
# Standard pipeline parameter
rdp=/mnt/bai/public/ref/rdp_gold.fa # rdp gold database, fro remove chimera
gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta # greengene bacterial 16S database
gg_seq=/mnt/bai/public/ref/gg_13_8_otus/rep_set/97_otus.fasta
gg_tax=/mnt/bai/public/ref/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
#gg_seq=/mnt/bai/public/ref/rdp/Bacteria_Archaea_seq.fa
#gg_tax=/mnt/bai/public/ref/rdp/Bacteria_Archaea_tax.txt
log=result/readme.log # log file for basic statistics
bc1=0 # forword barcode length, default 0
bc2=6 # forword barcode length, default 0
quality=19 # base quality, accurate > 99%; 29 means 99.9%
bt=6 # barcode type, usually length equal barcode 1 add barcode 2
primer5=AACMGGATTAGATACCCKG # 5` primer used for 16S
primer3=GGAAGGTGGGGATGACGT # 3` primer used for 16S, must reverse compliment
min_len=300 # min length, recommend 300 for bacterial 16S and 220 for ITS
thre_count=3000 # sample min count, filter samples less than thre_count
minuniquesize=2 # min count of unique reads
sim=0.97 # similarity of cluster OTU
p=32 # threads number used: 32
tax_per=0.005 # filter OTU percentage > 0.5% for draw taxonomy and phylogenetic tree, 0.1% about 150 OTU is too much to show
method=rdp # rdp, blast, rtax, mothur, uclust, sortmerna , default=uclust

# OTU taxonomy and abundance filter parameter
thre=0.001 # threshold of filter low abundance OTU
taxonomy=p__Cyanobacteria,p__Chloroflexi # filter some phylum
result=result_k1-c # result based on filter OTU table

align_seqs.py -i temp/otus_rdp.fa -t $gg_align -o temp/aligned/
fasta_subtraction.pl -i temp/otus_rdp.fa -d temp/aligned/otus_rdp_failures.fasta -o temp/otus_rdp_align.fa
echo 'Remove non-bac seq by align_seqs.py:' >> $log
grep '>' -c temp/otus_rdp_align.fa >> $log


# 8. Taxonomy assignment
assign_taxonomy.py -i result/rep_seqs.fa -r $gg_seq -t $gg_tax -m ${method} -o result
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt # format for R read
mv result/rep_seqs_tax_assignments.log temp/rep_seqs_tax_assignments.log
biom add-metadata -i result/otu_table.biom --observation-metadata-fp result/rep_seqs_tax_assignments.txt -o result/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy # add taxonomy to biom
biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy
summarize_taxa.py -i result/otu_table_tax.biom -o result/sum_taxa # summary each level percentage
rm result/sum_taxa/*.biom
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* # format for R read

# 9. Phylogeny tree
clustalo -i result/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=$p
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree # generate tree by FastTree

# 10. Alpha diversity
rarefaction=`head -n 7 result/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
single_rarefaction.py -i result/otu_table.biom -o temp/otu_table_rare.biom -d $rarefaction
alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 11. Beta diversity
normalize_table.py -i result/otu_table.biom -o temp/otu_table_css.biom -a CSS
biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
sed -i 's/^\t//g' result/beta/*

# 12. Taxonomy tree - GraPhlAn
filter_otus_from_otu_table.py --min_count_fraction $tax_per -i result/otu_table.biom -o temp/tax_otu_table.biom
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/tax_otu_table.biom 
echo "Number of OTU abundance > $tax_per :" >> $log
grep -c '>' temp/tax_rep_seqs.fa >> $log
grep '>' temp/tax_rep_seqs.fa|sed 's/>//g' > temp/tax_rep_seqs.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs.id|cut -f 2-3|grep 's__'|sed 's/; */\|/g' > temp/tax_full_anno.txt # rdp after ; no blank
echo "Number of OTU abundance > $tax_per with fully annotation :" >> $log
wc -l temp/tax_full_anno.txt >> $log
echo "Number of OTU abundance > $tax_per with fully annotation unique:" >> $log
sort temp/tax_full_anno.txt|cut -f 1|uniq|wc -l >> $log
format_taxonomy2lefse.pl -i temp/tax_full_anno.txt -o temp/tax_lefse.txt 
## order
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_order.tree --annotation temp/tax_order.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 4 --min_clade_size 1 --min_font_size 5
graphlan_annotate.py --annot temp/tax_order.annot temp/tax_order.tree temp/tax_order.xml
sed -i 's/ref="A:1">o  /ref="A:1">/g' temp/tax_order.xml
graphlan.py --dpi 300 temp/tax_order.xml result/tax_order.pdf --external_legends
mv result/tax_order_legend.pdf temp/ 
## family
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_family.tree --annotation temp/tax_family.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 5 --min_clade_size 1 --min_font_size 4
graphlan_annotate.py --annot temp/tax_family.annot temp/tax_family.tree temp/tax_family.xml
sed -i 's/ref="A:1">f  /ref="A:1">/g' temp/tax_family.xml
graphlan.py --dpi 300 temp/tax_family.xml result/tax_family.pdf --external_legends
mv result/tax_family_legend.pdf temp/ 
## genus
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_genus.tree --annotation temp/tax_genus.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 6 --min_clade_size 1 --min_font_size 3
graphlan_annotate.py --annot temp/tax_genus.annot temp/tax_genus.tree temp/tax_genus.xml
sed -i 's/ref="A:1">g  /ref="A:1">/g' temp/tax_genus.xml
graphlan.py --dpi 300 temp/tax_genus.xml result/tax_genus.pdf --external_legends
mv result/tax_genus_legend.pdf temp/ 

# 13. Phylogenetic tree - ggtree
clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=$p
make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/; /\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax

# R visualization
diversity.sh -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D 1 -o result -m FALSE -p TRUE # select group1 - 5, output, no merge group, pair compare
diversity.sh -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D 2 -o result -m FALSE -p TRUE # select group2 - 5
diversity.sh -A genotype -C batch -D 2 -o result -m TRUE # select group2 - all

taxonomy_lm.sh -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D 1 -o result -m FALSE -p TRUE

DAOTU_egr_p.sh -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D 1 -o result -m FALSE -p TRUE 

Rscript diversity.r # Draw alpha, beta and Constrain PCoA
Rscript taxonomy.R # Draw barplot+error bar, stack plot
Rscript DEOTU.R # Draw volcano, manhattan, heatmap, venn


# Optional 1. Filter taxonomy and low abundance OTU

# QIIME filter OTU table
mkdir $result
cd result
Rscript filter_otu_table.R
cd ..
filter_otus_from_otu_table.py -i result/otu_table_tax.biom -o temp/k1.biom --otu_ids_to_exclude_fp result/otu_id_k1.txt --negate_ids_to_exclude
echo 'Summary of otu_table_k1, one of sample OTU > 0.1%:' >> $log
biom summarize-table -i temp/k1.biom  >> $log
filter_taxa_from_otu_table.py -i temp/k1.biom -o $result/otu_table.biom -n $taxonomy
echo 'Summary of otu_table_k1 remove:'$taxonomy >> $log
biom summarize-table -i $result/otu_table.biom >> $log
biom summarize-table -i $result/otu_table.biom > $result/otu_table.sum
filter_fasta.py -f result/rep_seqs.fa -o $result/rep_seqs.fa -b $result/otu_table.biom
ln $result/otu_table.biom $result/otu_table_tax.biom
summarize_taxa.py -i $result/otu_table_tax.biom -o $result/sum_taxa
rm $result/sum_taxa/*.biom
sed -i '/# Const/d;s/#OTU ID.//g' $result/sum_taxa/*
biom convert -i $result/otu_table.biom -o $result/otu_table.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU ID.//' $result/otu_table.txt 
cut -f 1 $result/otu_table.txt | tail -n+2 > temp/k1_t.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax.txt temp/k1_t.id > $result/rep_seqs_tax.txt

# 9. Phylogeny tree
clustalo -i $result/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=$p
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o $result/rep_seqs.tree # generate tree by FastTree

# 10. Alpha diversity
rarefaction=`head -n 7 $result/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
single_rarefaction.py -i $result/otu_table.biom -o temp/otu_table_rare.biom -d $rarefaction
alpha_diversity.py -i temp/otu_table_rare.biom -o $result/alpha.txt -t $result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 11. Beta diversity
normalize_table.py -i $result/otu_table.biom -o temp/otu_table_css.biom -a CSS
biom convert -i temp/otu_table_css.biom -o $result/otu_table_css.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' $result/otu_table_css.txt
beta_diversity.py -i temp/otu_table_css.biom -o $result/beta/ -t $result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
sed -i 's/^\t//g' $result/beta/*

# R visualization
cd $result
Rscript diversity.R # Draw alpha, beta and Constrain PCoA
Rscript taxonomy_bar.R # Draw barplot+error bar, stack plot
Rscript DEOTU.R # Draw volcano, manhattan, heatmap, venn

## Addtional 1. Network by CoNet
# OTU network
# Example of D:\work\test\CoNet, using DM1 WT DO2 of batch1 as test





