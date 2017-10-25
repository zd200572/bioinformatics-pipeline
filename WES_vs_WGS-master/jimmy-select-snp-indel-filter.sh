module load java/1.8.0_91
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta 
GATK=/home/jianmingzeng/biosoft/GATK/GenomeAnalysisTK.jar 
TMPDIR=/home/jianmingzeng/tmp/software 
sample=$1
## for SNP 
: '
'
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants  -R $GENOME  \
-selectType SNP \
-V ${sample}_raw.snps.indels.vcf -o ${sample}_raw_snps.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantFiltration -R $GENOME  \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  \
--filterName "my_snp_filter" \
-V ${sample}_raw_snps.vcf  -o ${sample}_filtered_snps.vcf   
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants -R $GENOME  \
--excludeFiltered \
-V ${sample}_filtered_snps.vcf  -o  ${sample}_filtered_PASS.snps.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantEval -R $GENOME  \
-eval ${sample}_filtered_PASS.snps.vcf -o  ${sample}_filtered_PASS.snps.vcf.eval
## for  INDEL 
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants  -R $GENOME \
-selectType INDEL  \
-V ${sample}_raw.snps.indels.vcf   -o ${sample}_raw_indels.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantFiltration -R $GENOME  \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"  \
--filterName "my_indel_filter" \
-V ${sample}_raw_indels.vcf  -o ${sample}_filtered_indels.vcf   
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants -R $GENOME  \
--excludeFiltered \
-V ${sample}_filtered_indels.vcf  -o  ${sample}_filtered_PASS.indels.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantEval -R $GENOME  \
-eval ${sample}_filtered_PASS.indels.vcf -o  ${sample}_filtered_PASS.indels.vcf.eval
























GenomeAnalysisTK -R chr1.fa \
-T PrintReads \
-BQSR A.recal \
-o A.chr1.sort.dedup.mate.relaign.recal.bam \
-I A.chr1.sort.dedup.mate.relaign.bam

GenomeAnalysisTK -R chr1.fa -T UnifiedGenotyper -glm BOTH -D ../WES/ -metrics A.snps.metrics -stand_call_conf 50.0 -stand_emit_conf 10.0 -o A.ug.vcf 
-I A.chr1.sort.dedup.mate.relaign.recal.bam

GenomeAnalysisTK -T UnifiedGenotyper -R chr1.fa -I  A.chr1.sort.dedup.mate.bam -L 20 -glm BOTH -stand_call_conf 30 -stand_emit_conf 10 -o raw_ug_variants.vcf 