# ��������Ŀ¼������,-p����Ϊ����ļ��д��ڲ�����
mkdir -p example_PE250
cd example_PE250
# ����ʱ�ļ��ͽ����Ŀ¼
mkdir -p temp result
��
# �������������Ҫ��ʮ��
source activate qiime2-2017.7

# ����Ƿ�װ�ɹ�����������������ɹ�
qiime --help

# �رչ�������������ʱ�رգ���Ȼ������������ܻ����
source deactivate
��
# ��֤ʵ������Ƿ��д���
validate_mapping_file.py -m mappingfile.txt
#˫�����ݺϲ�Ϊ�����ļ�
join_paired_ends.py -f PE250_1.fq.gz -r PE250_2.fq.gz -m fastq-join -o temp/PE250_join
# ��ȡbarcode
extract_barcodes.py -f temp/PE250_join/fastqjoin.join.fastq -m mappingfile.txt -o temp/PE250_barcode -c barcode_paired_stitched --bc1_len 0 --bc2_len 6 -a --rev_comp_bc2

reads.fastq # �����ļ�����barcode��Ӧ���������η���
# �ʿؼ���Ʒ���
split_libraries_fastq.py -i temp/PE250_barcode/reads.fastq -b temp/PE250_barcode/barcodes.fastq  -m mappingfile.txt -o temp/PE250_split/ -q 20 --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type 6

#cutadapt�г�˫�����Ｐ���ȿ���
cutadapt -g AACMGGATTAGATACCCKG -a GGAAGGTGGGGATGACGT -e 0.15 -m 300 --discard-untrimmed temp/PE250_split/seqs.fna -o temp/PE250_P5.fa

# ��ʽת��
sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' temp/PE250_P5.fa > temp/seqs_usearch.fa

# ����ȥ����
./usearch10 -fastx_uniques temp/seqs_usearch.fa -fastaout temp/seqs_unique.fa -minuniquesize 2 -sizeout

# ����OTU
 ./usearch10.0.240_i86linux32 -cluster_otus temp/seqs_unique.fa -otu
 s temp/otus.fa -uparseout temp/uparse.txt -relabel Otu

# �鿴OTU����
grep '>' -c temp/otus.fa
######04:26 85Mb    100.0% 5493 OTUs, 9086 chimeras

# ����Usearch8�Ƽ��Ĳο����ݿ�RDP
wget http://drive5.com/uchime/rdp_gold.fa
# ����RDP���ݿ�ȶ�ȥ����֪���е�Ƕ����
./usearch10 -uchime2_ref temp/otus.fa -db rdp_gold.fa -chimeras temp/otus_chimeras.fa -notmatched temp/otus_rdp.fa -uchimeout temp/otus_rdp.uchime -strand plus -mode sensitive -threads 4

# ����Greengene�������ݿ⣬320MB
wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# ��ѹ���ݰ����С3.4G
tar xvzf gg_13_8_otus.tar.gz
# ��OTU��97%���ƾ���Ĵ��������ж����бȶԣ���Լ8min
time align_seqs.py -i temp/otus_non_chimera.fa -t gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
# �޷��ȶ�ϸ��������
grep -c '>' temp/aligned/otus_non_chimera_failures.fasta # 1860
# ��ò���ϸ����OTU ID
grep '>' temp/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_non_chimera_failures.id
# ���˷�ϸ������
filter_fasta.py -f temp/otus_non_chimera.fa -o temp/otus_rdp_align.fa -s temp/aligned/otus_non_chimera_failures.id -n
# ���������ڻ��ж���OTU:975
grep '>' -c temp/otus_rdp_align.fa

# ���Ƕ���������ID
grep '>' temp/otus_chimeras.fa|sed 's/>//g' > temp/otus_chimeras.id
# �޳�Ƕ���������
filter_fasta.py -f temp/otus.fa -o temp/otus_non_chimera.fa -s temp/otus_chimeras.id -n
# ����Ƿ�ΪԤ�ڵ���������
grep '>' -c temp/otus_non_chimera.fa # 2835

# ����Greengene�������ݿ⣬320MB
wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# ��ѹ���ݰ����С3.4G
tar xvzf gg_13_8_otus.tar.gz
# ��OTU��97%���ƾ���Ĵ��������ж����бȶԣ���Լ8min
time align_seqs.py -i temp/otus_non_chimera.fa -t gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
# �޷��ȶ�ϸ��������
grep -c '>' temp/aligned/otus_non_chimera_failures.fasta # 1860
# ��ò���ϸ����OTU ID
grep '>' temp/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_non_chimera_failures.id
# ���˷�ϸ������
filter_fasta.py -f temp/otus_non_chimera.fa -o temp/otus_rdp_align.fa -s temp/aligned/otus_non_chimera_failures.id -n
# ���������ڻ��ж���OTU:975
grep '>' -c temp/otus_rdp_align.fa

# ������OTU����������հ�Ĵ��������У���Reference(��ѡ������ϰ��)
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
# ����OTU��
./usearch10 -usearch_global temp/seqs_usearch.fa -db result/rep_seqs.fa -otutabout temp/otu_table.txt -biomout temp/otu_table.biom -strand plus -id 0.97 -threads 4
# �����Ϣ 01:20 141Mb   100.0% Searching seqs_usearch.fa, 32.3% matched
# Ĭ��10�̣߳���ʱ1��20�룬��32.3%������ƥ�䵽OTU�ϣ���30�̷߳�����ʱ3��04�룬�����߳�Խ��Խ�죬�ַ�����Ҳ�Ǻܷ�ʱ���

#�����OTU����less temp/otu_table.txt�鿴һ��
#13 # ����ע��
assign_taxonomy.py -i result/rep_seqs.fa \
    -r gg_13_8_otus/rep_set/97_otus.fasta \
    -t gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
    -m rdp -o result
#14. OTU��ͳ�ơ���ʽת���������Ϣ

# �ı�OTU��ת��ΪBIOM���������
biom convert -i temp/otu_table.txt \
    -o result/otu_table.biom \
    --table-type="OTU table" --to-json
# ���������Ϣ��OTU�����һ�У�����Ϊtaxonomy
biom add-metadata -i result/otu_table.biom \
    --observation-metadata-fp result/rep_seqs_tax_assignments.txt \
    -o result/otu_table_tax.biom \
    --sc-separated taxonomy --observation-header OTUID,taxonomy 
# ת��biomΪtxt��ʽ����������ע�ͣ�����ɶ�
biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy

# �鿴OTU��Ļ�����Ϣ����Ʒ��OUT����ͳ��
biom summarize-table -i result/otu_table_tax.biom -o result/otu_table_tax.sum
less result/otu_table_tax.sum#�鿴һ��
#15. OTU��ɸѡ

# ����Ʒ���������ˣ�ѡ��counts>3000����Ʒ
filter_samples_from_otu_table.py -i result/otu_table_tax.biom -o result/otu_table2.biom -n 3000
# �鿴���˺�����ֻ��25����Ʒ��975��OTU
biom summarize-table -i result/otu_table2.biom
# ��OTU��ȹ��ˣ�ѡ����Է�Ⱦ�ֵ�������֮һ��OTU
filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i result/otu_table2.biom -o result/otu_table3.biom
# �鿴���˺�����ֻ��25����Ʒ��346��OTU
biom summarize-table -i result/otu_table3.biom
# ת������biom��ʽOTU��Ϊ�ı�OTU���
biom convert -i result/otu_table3.biom -o result/otu_table4.txt --table-type="OTU table" --to-tsv
# OTU���ʽ��������R��ȡ
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt
# ɸѡ����OTU���ж�Ӧ��OTU����
filter_fasta.py -f result/rep_seqs.fa -b result/otu_table3.biom -o result/rep_seqs4.fa
#16. ����������

#�������ǻ��ڶ����бȶԵĽ������չʾ�ḻ����Ϣ�����ǽ���R��ͼ����ϸ������˴�ֻ�ǽ���������Alpha, Beta�����Է����������ļ���

# clustalo�����бȶԣ����û���밲װClustal Omega
clustalo -i result/rep_seqs4.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=4
# ɸѡ����б������кͱ�����
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/rep_seqs_align_pfiltered#only very short conserved region saved
# ����fasttree����
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fa/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree 
# generate tree by FastTree


#17. Alpha������

#Alpha�������Ǽ�����Ʒ��������ɣ����������ͷ����ά��Ϣ��������Ϳɼ�1����ͼ��Alpha�����ԣ��ϰ���Ҳ�������ҵ������Ķ�  

#Alpha�����Լ���ǰ��Ҫ��OTU����б�׼������Ϊ��ͬ������ȣ���⵽�����������᲻ͬ�����ǽ�OTU���س�������ͬ���������Թ�ƽ�Ƚϸ���Ʒ�������������������£�

# �鿴��Ʒ����������Сֵ
biom summarize-table -i result/otu_table3.biom
# ������Сֵ�����س�����׼��
single_rarefaction.py -i result/otu_table3.biom -o temp/otu_table_rare.biom -d 2797
# ���㳣�õ�����Alpha������ָ��
alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree



#18. Beta������

#Beta�������Ǽ������Ʒ�����ͬ��ͬ��OTU��Ҳ��Ҫ��׼���������س���������ʧ����Ϣ̫�࣬������ͳ�ơ��˲�����ѡ��CSS��׼��������

# CSS��׼��OTU��
normalize_table.py -i result/otu_table3.biom -o temp/otu_table_css.biom -a CSS
# ת����׼��OTU��Ϊ�ı������ں��ڻ�ͼ
biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
# ɾ����������Ϣ������R��ȡ
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
# ����Beta������
beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
# Beta�����Ծ����ļ���������R��ȡ
sed -i 's/^\t//g' result/beta/*


#19. �����ַ��༶��������

#OTU��������Ҫ��ע����Ϣ������ע����Ϣ��ͨ��������ע����Ϣ��Ϊ7�����𣺽硢�š��١�Ŀ���ơ������֡�������С�ļ��𣬺�OTU���Ƶ��в���ͬ��
#���ǳ��˿��ԱȽ���Ʒ�����OTUˮƽ�����⣬�������о���ͬ���Ƽ����ϵĲ��죬�����Ƿ������Щ��ͬ�ı仯���ɡ�

#����ע�͵ļ�����з�����ܣ�������Excel��R�������������Ǻ��鷳�Ĺ��̡���������ʹ��QIIME�Դ� �Ľű�summarize_taxa.py��

# ������š��١�Ŀ���ơ������������з�����ܣ���Ӧ�����L2-L6
summarize_taxa.py -i result/otu_table3.biom -o result/sum_taxa
 # summary each level percentage
# �޸�һ���ı���ͷ���ʺ�R��ȡ�ı���ʽ
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* 
# format for R read
# ����Ϊ���鿴���
less -S result/sum_taxa/otu_table3_L2.txt


#20. ɸѡ��չʾ�Ľ�����

#�����������п�������Ư���Ľ�����������OTUͨ���ɰ���ǧ�����ֱ��չʾ�Ǹ���������Ҳ�Ǽ���ġ�
#����̴��һЩͨ���ķ�����ɸѡ���ݣ������ó�Ư���Ľ�������

# ѡ��OTU���з�ȴ���0.1%��OTU
filter_otus_from_otu_table.py --min_count_fraction 0.001 -i result/otu_table3.biom -o temp/otu_table_k1.biom
# ��ö�Ӧ��fasta����
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/otu_table_k1.biom 
# ͳ������������104����һ��100�����Ҽ��д����ݵ�B�����ܶ����͸�����ɺ�ϸ��
grep -c '>' temp/tax_rep_seqs.fa # 104
# �����бȶ�
clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=4
# ����
make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
# ��ʽת��ΪR ggtree���õ���
sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
# �������ID
grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
# �����Щ���е�����ע�ͣ�����������ɫ��ʾ��ͬ������Ϣ
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/; /\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax


#21. ����

#��������һЩ�򵥵ĸ�ʽת����Ϊ����ͳ�Ʒ�����׼���ļ����Ժ�������ʲô���⣬�����ﲹ�䣬��CSDN�������ҿ����޸ģ����Ѵ�ҷ���������ʹ𰸷����ⲿ�֡�

# ��mappingfileת��ΪR�ɶ���ʵ�����
sed 's/#//' mappingfile.txt > result/design.txt
# ת���ı�otu_table��ʽΪR�ɶ�
sed '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table4.txt > result/otu_table.txt
# ת������ע����ϢΪ�Ʊ���ָ�������R��ȡ
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt








