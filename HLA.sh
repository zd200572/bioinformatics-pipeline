###########################
#HLA NGS analysis pipeline
###########################
#测序文件质控
ls *.fastq | while read id; do fastqc -t 4 $id; done 
#参考序列下载，索引制作
#mkdir reference && cd reference
'''
#下载
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
#解压
tar -zvf hg38.chromFa.tar.gz
#合并
cat *.fa > hg38.fa
#删除多余文件
rm chr*
'''
#bwa比对，先是
#cd ..
#bwa index reference/hg38.fa #建立索引
ls *gz | while read id; do bwa mem reference/hg38.fa $id > $id.sam; done #然后比对
ls *sam | while read id; do samtools view -S $id -b > $id.bam; done #格式转换
ls *sam.bam | while read id; do samtools sort $id $id.sorted.bam; done # 排序
ls *sorted.bam | while read id; do samtools index $id; done  # 索引
#xHLA 算法
ls *sorted.bam | while read id; do docker run -v `pwd`:`pwd` -w `pwd` humanlongevity/hla --sample_id test --input_bam_path $id --output_path test; done
################################
#end of pipeline
################################