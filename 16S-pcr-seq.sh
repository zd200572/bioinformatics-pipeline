# 建立工作目录并进入,-p参数为如果文件夹存在不报错
mkdir -p example_PE250
cd example_PE250
# 建临时文件和结果子目录
mkdir -p temp result
‘
# 激活工作环境，需要几十秒
source activate qiime2-2017.7

# 检查是否安装成功，弹出程序帮助即成功
qiime --help

# 关闭工作环境：不用时关闭，不然你其它程序可能会出错
source deactivate
’
# 验证实验设计是否有错误
validate_mapping_file.py -m mappingfile.txt
#双端数据合并为单个文件
join_paired_ends.py -f PE250_1.fq.gz -r PE250_2.fq.gz -m fastq-join -o temp/PE250_join
# 提取barcode
extract_barcodes.py -f temp/PE250_join/fastqjoin.join.fastq -m mappingfile.txt -o temp/PE250_barcode -c barcode_paired_stitched --bc1_len 0 --bc2_len 6 -a --rev_comp_bc2

reads.fastq # 序列文件，与barcode对应，用于下游分析
# 质控及样品拆分
split_libraries_fastq.py -i temp/PE250_barcode/reads.fastq -b temp/PE250_barcode/barcodes.fastq  -m mappingfile.txt -o temp/PE250_split/ -q 20 --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type 6

#cutadapt切除双端引物及长度控制
cutadapt -g AACMGGATTAGATACCCKG -a GGAAGGTGGGGATGACGT -e 0.15 -m 300 --discard-untrimmed temp/PE250_split/seqs.fna -o temp/PE250_P5.fa

# 格式转换
sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' temp/PE250_P5.fa > temp/seqs_usearch.fa

# 序列去冗余
./usearch10 -fastx_uniques temp/seqs_usearch.fa -fastaout temp/seqs_unique.fa -minuniquesize 2 -sizeout

# 聚类OTU
 ./usearch10.0.240_i86linux32 -cluster_otus temp/seqs_unique.fa -otu
 s temp/otus.fa -uparseout temp/uparse.txt -relabel Otu

# 查看OTU数量
grep '>' -c temp/otus.fa
######04:26 85Mb    100.0% 5493 OTUs, 9086 chimeras

# 下载Usearch8推荐的参考数据库RDP
wget http://drive5.com/uchime/rdp_gold.fa
# 基于RDP数据库比对去除已知序列的嵌合体
./usearch10 -uchime2_ref temp/otus.fa -db rdp_gold.fa -chimeras temp/otus_chimeras.fa -notmatched temp/otus_rdp.fa -uchimeout temp/otus_rdp.uchime -strand plus -mode sensitive -threads 4

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

# 获得嵌合体的序列ID
grep '>' temp/otus_chimeras.fa|sed 's/>//g' > temp/otus_chimeras.id
# 剔除嵌合体的序列
filter_fasta.py -f temp/otus.fa -o temp/otus_non_chimera.fa -s temp/otus_chimeras.id -n
# 检查是否为预期的序列数量
grep '>' -c temp/otus_non_chimera.fa # 2835

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

# 重命名OTU，这就是最终版的代表性序列，即Reference(可选，个人习惯)
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
# 生成OTU表
./usearch10 -usearch_global temp/seqs_usearch.fa -db result/rep_seqs.fa -otutabout temp/otu_table.txt -biomout temp/otu_table.biom -strand plus -id 0.97 -threads 4
# 结果信息 01:20 141Mb   100.0% Searching seqs_usearch.fa, 32.3% matched
# 默认10线程，用时1分20秒，有32.3%的序列匹配到OTU上；用30线程反而用时3分04秒，不是线程越多越快，分发任务也是很费时间的

#获得了OTU表，用less temp/otu_table.txt查看一下
#13 # 物种注释
assign_taxonomy.py -i result/rep_seqs.fa \
    -r gg_13_8_otus/rep_set/97_otus.fasta \
    -t gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
    -m rdp -o result
#14. OTU表统计、格式转换、添加信息

# 文本OTU表转换为BIOM：方便操作
biom convert -i temp/otu_table.txt \
    -o result/otu_table.biom \
    --table-type="OTU table" --to-json
# 添加物种信息至OTU表最后一列，命名为taxonomy
biom add-metadata -i result/otu_table.biom \
    --observation-metadata-fp result/rep_seqs_tax_assignments.txt \
    -o result/otu_table_tax.biom \
    --sc-separated taxonomy --observation-header OTUID,taxonomy 
# 转换biom为txt格式，带有物种注释：人类可读
biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy

# 查看OTU表的基本信息：样品，OUT数量统计
biom summarize-table -i result/otu_table_tax.biom -o result/otu_table_tax.sum
less result/otu_table_tax.sum#查看一下
#15. OTU表筛选

# 按样品数据量过滤：选择counts>3000的样品
filter_samples_from_otu_table.py -i result/otu_table_tax.biom -o result/otu_table2.biom -n 3000
# 查看过滤后结果：只有25个样品，975个OTU
biom summarize-table -i result/otu_table2.biom
# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i result/otu_table2.biom -o result/otu_table3.biom
# 查看过滤后结果：只有25个样品，346个OTU
biom summarize-table -i result/otu_table3.biom
# 转换最终biom格式OTU表为文本OTU表格
biom convert -i result/otu_table3.biom -o result/otu_table4.txt --table-type="OTU table" --to-tsv
# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt
# 筛选最终OTU表中对应的OTU序列
filter_fasta.py -f result/rep_seqs.fa -b result/otu_table3.biom -o result/rep_seqs4.fa
#16. 进化树构建

#进化树是基于多序列比对的结果，可展示丰富的信息，我们将在R绘图中详细解读。此处只是建树，用于Alpha, Beta多样性分析的输入文件。

# clustalo多序列比对，如果没有请安装Clustal Omega
clustalo -i result/rep_seqs4.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=4
# 筛选结果中保守序列和保守区
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/rep_seqs_align_pfiltered#only very short conserved region saved
# 基于fasttree建树
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fa/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree 
# generate tree by FastTree


#17. Alpha多样性

#Alpha多样性是计算样品内物种组成，包括数量和丰度两维信息。具体解释可见1箱线图：Alpha多样性，老板再也不操心我的文献阅读  

#Alpha多样性计算前需要对OTU表进行标准化，因为不同测序深度，检测到的物种数量会不同。我们将OTU表重抽样至相同数据量，以公平比较各样品的物种数量。方法如下：

# 查看样品的数据量最小值
biom summarize-table -i result/otu_table3.biom
# 基于最小值进行重抽样标准化
single_rarefaction.py -i result/otu_table3.biom -o temp/otu_table_rare.biom -d 2797
# 计算常用的四种Alpha多样性指数
alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree



#18. Beta多样性

#Beta多样性是计算各样品间的相同或不同，OTU表也需要标准化。采用重抽样方法丢失的信息太多，不利于统计。此步我们选择CSS标准化方法。

# CSS标准化OTU表
normalize_table.py -i result/otu_table3.biom -o temp/otu_table_css.biom -a CSS
# 转换标准化OTU表为文本，用于后期绘图
biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
# 删除表格多余信息，方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
# 计算Beta多样性
beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
# Beta多样性距离文件整理，方便R读取
sed -i 's/^\t//g' result/beta/*


#19. 按物种分类级别分类汇总

#OTU表中最重要的注释信息是物种注释信息。通常的物种注释信息分为7个级别：界、门、纲、目、科、属、种。种是最小的级别，和OTU类似但有不相同。
#我们除了可以比较样品和组间OTU水平差异外，还可以研究不同类似级别上的差异，它们是否存在那些共同的变化规律。

#按照注释的级别进行分类汇总，无论是Excel还R操作起来，都是很麻烦的过程。这里我们使用QIIME自带 的脚本summarize_taxa.py。

# 结果按门、纲、目、科、属五个级别进行分类汇总，对应结果的L2-L6
summarize_taxa.py -i result/otu_table3.biom -o result/sum_taxa
 # summary each level percentage
# 修改一下文本表头，适合R读取的表格格式
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* 
# format for R read
# 以门为例查看结果
less -S result/sum_taxa/otu_table3_L2.txt


#20. 筛选可展示的进化树

#我们在文章中看到几种漂亮的进化树，但是OTU通常成百上千，如果直接展示是根本看不清也是极丑的。
#下面教大家一些通常的方江来筛选数据，用于用成漂亮的进化树。

# 选择OTU表中丰度大于0.1%的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.001 -i result/otu_table3.biom -o temp/otu_table_k1.biom
# 获得对应的fasta序列
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/otu_table_k1.biom 
# 统计序列数量，104条，一般100条左右即有大数据的B格，又能读懂和更清规律和细节
grep -c '>' temp/tax_rep_seqs.fa # 104
# 多序列比对
clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=4
# 建树
make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
# 格式转换为R ggtree可用的树
sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
# 获得序列ID
grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
# 获得这些序列的物种注释，用于树上着色显示不同分类信息
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/; /\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax


#21. 其它

#其它都是一些简单的格式转换，为后其统计分析而准备文件。以后再碰到什么问题，在这里补充，在CSDN的文章我可以修改，并把大家反馈的问题和答案放在这部分。

# 将mappingfile转换为R可读的实验设计
sed 's/#//' mappingfile.txt > result/design.txt
# 转换文本otu_table格式为R可读
sed '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt > result/otu_table.txt
# 转换物种注释信息为制表符分隔，方便R读取
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt








