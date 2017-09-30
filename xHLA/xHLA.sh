###########################
#HLA NGS analysis pipeline
###########################
#测序文件质控

ls *.fastq | while read id; do fastqc -t 4 $id; done
ls *fastq | while read id; do bwa mem ../../reference_hg38_no_alt/chr6.fa $id > $id.sam; done #然后比对
ls *sam | while read id; do samtools view -S $id -b > $id.bam; done #格式转换
ls *sam.bam | while read id; do samtools sort $id -o $id.sorted.bam; done # 排序
ls *sorted.bam | while read id; do samtools index $id; done  # 索引

ls *sorted.bam | while read id; do sudo docker run -v `pwd`:`pwd` -w `pwd` humanlongevity/hla --sample_id test --input_bam_path $id --output_path test; done
