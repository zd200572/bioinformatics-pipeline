#cd fastq-join
#拼接
join_paired_ends.py -f /mnt/75/16s/$1 -r /mnt/75/16s/$2 -o fastq-join
#质控
split_libraries_fastq.py -q 19 --barcode_type not-barcoded --store_demultiplexed_fastq  -o temp --sample_ids test_i -i fastq-join/fastqjoin.join.fastq
## 格式转换
sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' temp/seqs.fna > temp/seqs_usearch.fa
# 序列去冗余
~/16s-sygene/usearch10 -fastx_uniques temp/seqs_usearch.fa -fastaout temp/seqs_unique.fa -minuniquesize 2 -sizeout
# 聚类OTU
~/16s-sygene/usearch10 -cluster_otus temp/seqs_unique.fa -otus temp/otus.fa -uparseout temp/uparse.txt -relabel Otu
# 查看OTU数量
grep '>' -c temp/otus.fa
# 下载Usearch8推荐的参考数据库RDP
#wget http://drive5.com/uchime/rdp_gold.fa
# 基于RDP数据库比对去除已知序列的嵌合体(可选)
#./usearch10 -uchime2_ref temp/otus.fa \
#    -db rdp_gold.fa \
#    -chimeras temp/otus_chimeras.fa \
#    -notmatched temp/otus_rdp.fa \
#    -uchimeout temp/otus_rdp.uchime \
#    -strand plus -mode sensitive -threads 4
# 获得嵌合体的序列ID
#grep '>' temp/otus_chimeras.fa|sed 's/>//g' > temp/otus_chimeras.id
# 剔除嵌合体的序列
#filter_fasta.py -f temp/otus.fa -o temp/otus_non_chimera.fa -s temp/otus_chimeras.id -n
# 检查是否为预期的序列数量
#grep '>' -c temp/otus_non_chimera.fa # 2835
#去除非细菌序列(可选)
# 下载Greengene最新数据库，320MB
#wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# 解压数据包后大小3.4G
#tar xvzf gg_13_8_otus.tar.gz
# 将OTU与97%相似聚类的代表性序列多序列比对，大约8min
#time align_seqs.py -i temp/otus.fa -t gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
# 无法比对细菌的数量
#grep -c '>' temp/aligned/otus_non_chimera_failures.fasta # 1860
# 获得不像细菌的OTU ID
#grep '>' temp/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_non_chimera_failures.id
# 过滤非细菌序列
#filter_fasta.py -f temp/otus_non_chimera.fa -o temp/otus_rdp_align.fa -s temp/aligned/otus_non_chimera_failures.id -n
# 看我们现在还有多少OTU:975
#grep '>' -c temp/otus_rdp_align.fa
#####产生代表性序列和OTU表
# 重命名OTU，这就是最终版的代表性序列，即Reference(可选，个人习惯)
#awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
# 生成OTU表
~/16s-sygene/usearch10 -usearch_global temp/seqs_usearch.fa -db temp/otus.fa -otutabout temp/otu_table.txt -biomout temp/otu_table.biom -strand plus -id 0.97 -threads 10
# 结果信息 01:20 141Mb   100.0% Searching seqs_usearch.fa, 32.3% matched
# 默认10线程，用时1分20秒，有32.3%的序列匹配到OTU上；用30线程反而用时3分04秒，不是线程越多越快，分发任务也是很费时间的
# 物种注释
assign_taxonomy.py -i temp/otus.fa \
    -r /home/biolinux/16s-sygene/gg_13_8_otus/rep_set/97_otus.fasta \
    -t /home/biolinux/16s-sygene/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
    -m rdp -o result
# 文本OTU表转换为BIOM：方便操作
biom convert -i temp/otu_table.txt \
    -o result/otu_table.biom \
    --table-type="OTU table" --to-json
# 添加物种信息至OTU表最后一列，命名为taxonomy
biom add-metadata -i result/otu_table.biom \
    --observation-metadata-fp result/otus_tax_assignments.txt \
    -o result/otu_table_tax.biom \
    --sc-separated taxonomy --observation-header OTUID,taxonomy 
# 转换biom为txt格式，带有物种注释：人类可读
biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy
# 查看OTU表的基本信息：样品，OUT数量统计
biom summarize-table -i result/otu_table_tax.biom -o result/otu_table_tax.sum
# 按样品数据量过滤：选择counts>3000的样品
filter_samples_from_otu_table.py -i result/otu_table_tax.biom -o result/otu_table2.biom -n 3000
# 查看过滤后结果：只有25个样品，975个OTU
biom summarize-table -i result/otu_table2.biom
# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i result/otu_table2.biom -o result/otu_table3.biom
# 查看过滤后结果：只有25个样品，346个OTU
biom summarize-table -i result/otu_table3.biom
# 转换最终biom格式OTU表为文本OTU表格
biom convert -i result/otu_table3.biom -o result/otu_table3.txt --table-type="OTU table" --to-tsv
# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table3.txt
# 筛选最终OTU表中对应的OTU序列
filter_fasta.py -f temp/otus.fa -b result/otu_table3.biom -o result/otus3.fa
# 结果按门、纲、目、科、属五个级别进行分类汇总，对应结果的L2-L6
summarize_taxa.py -i result/otu_table3.biom -o result/sum_taxa # summary each level percentage
# 修改一下文本表头，适合R读取的表格格式
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* # format for R read
# 以门为例查看结果
#less -S result/sum_taxa/otu_table3_L2.txt
cd ..