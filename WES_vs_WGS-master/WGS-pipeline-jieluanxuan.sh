#WGS-pipeline-jieluanxuan.sh
bwa mem -t 4 -R '@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:sample_name' /path/to/human.fasta read_1.fq.gz read_2.fq.gz | samtools view -S -b - > sample_name.bam

 time samtools sort -@ 4 -m 4G -O bam -o sample_name.sorted.bam sample_name.bam

 java -jar picard.jar MarkDuplicates \ 
  I=sample_name.sorted.bam \
  O=sample_name.sorted.markdup.bam \
  M=sample_name.markdup_metrics.txt

  #samtools merge <out.bam> <in1.bam> [<in2.bam> ... <inN.bam>]
  #多样本可选

  ava -jar /path/to/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.bam \
 -known /path/to/gatk/bundle/1000G_phase1.indels.b37.vcf \
 -known /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
 -o sample_name.IndelRealigner.intervals

java -jar /path/to/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.bam \
 -known /path/to/gatk/bundle/1000G_phase1.indels.b37.vcf \
 -known /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
 -o sample_name.sorted.markdup.realign.bam \
 --targetIntervals sample_name.IndelRealigner.intervals
#重新校正碱基质量值（BQSR）
 java -jar /path/to/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R /path/to/human.fasta \
    -I sample_name.sorted.markdup.realign.bam \
    —knownSites /path/to/gatk/bundle/1000G_phase1.indels.b37.vcf \
    —knownSites /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
    —knownSites /path/to/gatk/bundle/dbsnp_138.b37.vcf \
    -o sample_name.recal_data.table

java -jar /path/to/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R /path/to/human.fasta \
    -I sample_name.sorted.markdup.realign.bam \
    —BQSR sample_name.recal_data.table \
    -o sample_name.sorted.markdup.realign.BQSR.bam

#变异检测--单样本 -l 可以添加多样本，耗时 -L可以指定染色体
java -jar /path/to/GenomeAnalysisTK.jar \
 -T HaplotypeCaller \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.realign.BQSR.bam \
 -D /path/to/gatk/bundle/dbsnp_138.b37.vcf \
 -stand_call_conf 50 \ 
 -A QualByDepth \ 
 -A RMSMappingQuality \ 
 -A MappingQualityRankSumTest \ 
 -A ReadPosRankSumTest \ 
 -A FisherStrand \ 
 -A StrandOddsRatio \ 
 -A Coverage
 -o sample_name.HC.vcf
 #多样本--gvcf
 java -jar /path/to/GenomeAnalysisTK.jar \
 -T HaplotypeCaller \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.realign.BQSR.bam \
 --emitRefConfidence GVCF \
 -o sample_name.g.vcf
 #调用GenotypeGVCFs完成变异calling
 java -jar /path/to/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs \
 -R /path/to/human.fasta \
 --variant sample_name.g.vcf \
 -o sample_name.HC.vcf

 #变异检测质控和过滤（VQSR）
 ## SNP Recalibratorjava -jar /path/to/GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R reference.fasta \
   -input sample_name.HC.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /path/to/gatk/bundle/hapmap_3.3.b37.vcf \ 
   -resource:omini,known=false,training=true,truth=false,prior=12.0 /path/to/gatk/bundle/1000G_omni2.5.b37.vcf \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 /path/to/gatk/bundle/1000G_phase1.snps.high_confidence.b37.vcf \ 
   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 /path/to/gatk/bundle/dbsnp_138.b37.vcf \ 
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \ 
   -mode SNP \ 
   -recalFile sample_name.HC.snps.recal \
   -tranchesFile sample_name.HC.snps.tranches \ 
   -rscriptFile sample_name.HC.snps.plots.R

java -jar /path/to/GenomeAnalysisTK.jar -T ApplyRecalibration \
   -R  human_g1k_v37.fasta \
   -input sample_name.HC.vcf \ 
   --ts_filter_level 99.5 \ 
   -tranchesFile sample_name.HC.snps.tranches \ 
   -recalFile sample_name.HC.snps.recal \
   -mode SNP \
   -o sample_name.HC.snps.VQSR.vcf## Indel Recalibratorjava -jar /path/to/GenomeAnalysisTK.jar -T VariantRecalibrator \
   -R  human_g1k_v37.fasta \
   -input sample_name.HC.snps.VQSR.vcf \
   -resource:mills,known=true,training=true,truth=true,prior=12.0 /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
   -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -recalFile sample_name.HC.snps.indels.recal \
   -tranchesFile sample_name.HC.snps.indels.tranches \
   -rscriptFile sample_name.HC.snps.indels.plots.R

java -jar /path/to/GenomeAnalysisTK.jar -T ApplyRecalibration \ 
   -R human_g1k_v37.fasta\
   -input sample_name.HC.snps.VQSR.vcf \
   --ts_filter_level 99.0 \
   -tranchesFile sample_name.HC.snps.indels.tranches \
   -recalFile sample_name.HC.snps.indels.recal \
   -mode INDEL \
   -o sample_name.HC.snps.indels.VQSR.vcf

