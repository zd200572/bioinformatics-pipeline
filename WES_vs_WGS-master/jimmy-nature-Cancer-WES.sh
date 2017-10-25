#jimmy-nature-Cancer-WES.sh
#需要安装的软件
conda install -c bioconda bedtools
conda install -c bioconda bwa
conda install -c bioconda samtools
cd ~/biosoft
mkdir sratoolkit &&  cd sratoolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar zxvf sratoolkit.current-centos_linux64.tar.gz
~/biosoft/sratoolkit/sratoolkit.2.8.2-1-centos_linux64/bin/fastdump -h
## https://sourceforge.net/projects/picard/
## https://github.com/broadinstitute/picard
cd ~/biosoft
mkdir picardtools &&  cd picardtools
wget http://ncu.dl.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip
unzip picard-tools-1.119.zip
mkdir 2.9.2 && cd 2.9.2 
wget https://github.com/broadinstitute/picard/releases/download/2.9.2/picard.jar
cd ~/biosoft
## https://sourceforge.net/projects/varscan/files/
## http://varscan.sourceforge.net/index.html
mkdir VarScan  &&  cd VarScan  
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar 
cd ~/biosoft
mkdir SnpEff &&  cd SnpEff
##    http://snpeff.sourceforge.net/
##    http://snpeff.sourceforge.net/SnpSift.html 
##    http://snpeff.sourceforge.net/SnpEff_manual.html
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip 
## java -jar snpEff.jar download GRCh37.75
## java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.75 example.vcf > example_snpeff.vcf
unzip snpEff_latest_core.zip
#下载
cat npc.sra.txt | cut -f 4|while read id
do echo $id
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP035/SRP035573/$id/$id.sra
done 
#转换
cat npc.sra.txt | while read id
do
array=($id)
echo  ${array[3]}.sra  ${array[7]} 
~/biosoft/sratoolkit/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --gzip --split-3 -A \  
${array[7]}  ${array[3]}.sra 
done 
#质控
ls *.gz |xargs ~/biosoft/fastqc/FastQC/fastqc -o ./ -t 5 
#WES的标准SNP-calling流程
module load java/1.8.0_91
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
INDEX=/home/jianmingzeng/reference/index/bwa/human_g1k_v37
GATK=/home/jianmingzeng/biosoft/GATK/GenomeAnalysisTK.jar
PICARD=/home/jianmingzeng/biosoft/picardtools/2.9.2/picard.jar
DBSNP=/home/jianmingzeng/annotation/variation/human/dbSNP/All_20160601.vcf.gz
SNP=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
INDEL=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
KG_SNP=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
Mills_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
KG_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.indels.b37.vcf
TMPDIR=/home/jianmingzeng/tmp/software
## samtools and bwa are in the environment
## samtools Version: 1.3.1 (using htslib 1.3.1)
## bwa Version: 0.7.15-r1140
#arr=($1)
#fq1=${arr[0]}
#fq2=${arr[1]}
#sample=${arr[2]}
fq1=$1
fq2=$2
sample=$3
#####################################################
################ Step 1 : Alignment #################
#####################################################
start=$(date +%s.%N)
echo bwa `date`
bwa mem -t 5 -M  -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina" $INDEX $fq1 $fq2 1>$sample.sam 2>/dev/null 
echo bwa `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for BWA : %.6f seconds" $dur
echo 
#####################################################
################ Step 2: Sort and Index #############
#####################################################
start=$(date +%s.%N)
echo SortSam `date`
java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$sample.sam OUTPUT=$sample.bam
samtools index $sample.bam
echo SortSam `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for SortSam : %.6f seconds" $dur
echo 
rm $sample.sam 
#####################################################
################ Step 3: Basic Statistics ###########
#####################################################
start=$(date +%s.%N)
echo stats `date`
samtools flagstat $sample.bam > ${sample}.alignment.flagstat
samtools stats  $sample.bam > ${sample}.alignment.stat
echo plot-bamstats -p ${sample}_QC  ${sample}.alignment.stat
echo stats `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for Basic Statistics : %.6f seconds" $dur
echo 
#####################################################
####### Step 4: multiple filtering for bam files ####
#####################################################
###MarkDuplicates###
start=$(date +%s.%N)
echo MarkDuplicates `date`
java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD MarkDuplicates \
INPUT=$sample.bam OUTPUT=${sample}_marked.bam METRICS_FILE=$sample.metrics
echo MarkDuplicates `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for MarkDuplicates : %.6f seconds" $dur
echo 
rm $sample.bam  
###FixMateInfo###
start=$(date +%s.%N)
echo FixMateInfo `date`
java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD FixMateInformation \
INPUT=${sample}_marked.bam OUTPUT=${sample}_marked_fixed.bam SO=coordinate
samtools index ${sample}_marked_fixed.bam
echo FixMateInfo `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for FixMateInfo  : %.6f seconds" $dur
echo 
rm ${sample}_marked.bam 
#####################################################
####### Step 5: gatk process bam files ##############
#####################################################
### SplitNCigar ###
start=$(date +%s.%N)
echo SplitNCigar `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SplitNCigarReads \
 -R $GENOME  -I ${sample}_marked_fixed.bam  -o ${sample}_marked_fixed_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 
#--fix_misencoded_quality_scores
## --fix_misencoded_quality_scores only if phred 64 
echo SplitNCigar `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for SplitNCigar : %.6f seconds" $dur
echo 
rm ${sample}_marked_fixed.bam 
# rm ${sample}.sam ${sample}.bam ${sample}_marked.bam ${sample}_marked_fixed.bam
###RealignerTargetCreator###
start=$(date +%s.%N)
echo RealignerTargetCreator `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T RealignerTargetCreator \
-I ${sample}_marked_fixed_split.bam -R $GENOME -o ${sample}_target.intervals \
-known $Mills_indels -known $KG_indels -nt 5
echo RealignerTargetCreator `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for RealignerTargetCreator : %.6f seconds" $dur
echo 
###IndelRealigner###
start=$(date +%s.%N)
echo IndelRealigner `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T IndelRealigner \
-I ${sample}_marked_fixed_split.bam  -R $GENOME -targetIntervals ${sample}_target.intervals \
-o ${sample}_realigned.bam -known $Mills_indels -known $KG_indels
echo IndelRealigner `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for IndelRealigner : %.6f seconds" $dur
echo 
rm ${sample}_marked_fixed_split.bam
###BaseRecalibrator###
start=$(date +%s.%N)
echo BaseRecalibrator `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T BaseRecalibrator \
-I ${sample}_realigned.bam -R $GENOME -o ${sample}_temp.table -knownSites $DBSNP
echo BaseRecalibrator `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for BaseRecalibrator : %.6f seconds" $dur
echo 
###PrintReads###
start=$(date +%s.%N)
echo PrintReads `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T PrintReads \
-R $GENOME -I ${sample}_realigned.bam -o ${sample}_recal.bam -BQSR ${sample}_temp.table
samtools index ${sample}_recal.bam
echo PrintReads `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for PrintReads : %.6f seconds" $dur
echo 
rm  ${sample}_realigned.bam
chmod uga=r   ${sample}_recal.bam 
#####################################################
############## Step 6: gatk call snp/indel##########
#####################################################
### 
start=$(date +%s.%N)
echo HaplotypeCaller `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T HaplotypeCaller  \
-R $GENOME -I ${sample}_recal.bam --dbsnp $DBSNP  \
-stand_emit_conf 10 -o  ${sample}_raw.snps.indels.vcf
echo HaplotypeCaller `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for HaplotypeCaller : %.6f seconds" $dur
echo 
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants  -R $GENOME  \
-selectType SNP \
-V ${sample}_raw.snps.indels.vcf -o ${sample}_raw_snps.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants  -R $GENOME \
-selectType INDEL  \
-V ${sample}_raw.snps.indels.vcf   -o ${sample}_raw_indels.vcf
## 
:'
'
## for SNP
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
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantFiltration -R $GENOME  \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"  \
--filterName "my_indel_filter" \
-V ${sample}_raw_indels.vcf  -o ${sample}_filtered_indels.vcf   
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T SelectVariants -R $GENOME  \
--excludeFiltered \
-V ${sample}_filtered_indels.vcf  -o  ${sample}_filtered_PASS.indels.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T VariantEval -R $GENOME  \
-eval ${sample}_filtered_PASS.indels.vcf -o  ${sample}_filtered_PASS.indels.vcf.eval
#somatic mutation calling流程
reference=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
TMPDIR=/home/jianmingzeng/tmp/software
normal_bam=NPC10F-N_recal.bam
tumor_bam=NPC10F-T_recal.bam
sample=NPC10F
#####################################################
################### Step : Run VarScan #############
#####################################################
normal_pileup="samtools mpileup -q 1 -f $reference $normal_bam";
tumor_pileup="samtools mpileup -q 1 -f $reference $tumor_bam";
# Next, issue a system call that pipes input from these commands into VarScan :
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g  -jar ~/biosoft/VarScan/VarScan.v2.3.9.jar \
somatic <($normal_pileup) <($tumor_pileup) $sample
java -jar ~/biosoft/VarScan/VarScan.v2.3.9.jar processSomatic $sample.snp
#####################################################
################### Step : Run Mutect2 #############
##################################################### 
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK  -T MuTect2 \
-R $GENOME -I:tumor $tumor_bam  -I:normal $normal_bam \
--dbsnp  $DBSNP   -o ${sample}-mutect2.vcf
#####################################################
################### Step : Run Muse#################
#####################################################
~/biosoft/muse/muse call -O $sample -f $reference $tumor_bam $normal_bam
~/biosoft/muse/muse sump -I ${sample}.MuSE.txt -E –O ${sample}.vcf –D $DBSNP
