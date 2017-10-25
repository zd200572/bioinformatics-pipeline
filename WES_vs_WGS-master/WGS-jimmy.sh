#WGS-jimmy.sh
#软件安装
## Download and install BWA
cd ~/biosoft
mkdir bwa &&  cd bwa
#http://sourceforge.net/projects/bio-bwa/files/
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2 
tar xvfj bwa-0.7.15.tar.bz2 # x extracts, v is verbose (details of what it is doing), f skips prompting for each individual file, and j tells it to unzip .bz2 files
cd bwa-0.7.15
make
## Download and install samtools
## http://samtools.sourceforge.net/
## http://www.htslib.org/doc/samtools.html
cd ~/biosoft
mkdir samtools &&  cd samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 
tar xvfj samtools-1.3.1.tar.bz2 
cd samtools-1.3.1 
./configure --prefix=/home/jianmingzeng/biosoft/myBin
make 
make install 
~/biosoft/myBin/bin/samtools --help
~/biosoft/myBin/bin/plot-bamstats --help
cd htslib-1.3.1
./configure --prefix=/home/jianmingzeng/biosoft/myBin
make 
make install
~/biosoft/myBin/bin/tabix 
## Download and install picardtools
## https://sourceforge.net/projects/picard/
## https://github.com/broadinstitute/picard
cd ~/biosoft
mkdir picardtools &&  cd picardtools
wget http://ncu.dl.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip
unzip picard-tools-1.119.zip
mkdir 2.9.2 && cd 2.9.2 
wget https://github.com/broadinstitute/picard/releases/download/2.9.2/picard.jar
## GATK 需要自行申请下载，不能公开

#必备数据
cd ~/reference
mkdir -p  genome/human_g1k_v37  && cd genome/human_g1k_v37 
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/ 
nohup wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  &
gunzip human_g1k_v37.fasta.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/README.human_g1k_v37.fasta.txt
java -jar ~/biosoft/picardtools/picard-tools-1.119/CreateSequenceDictionary.jar R=human_g1k_v37.fasta O=human_g1k_v37.dict
cd ~/reference
mkdir -p index/bwa && cd index/bwa   ~/reference/index/bwa/human_g1k_v37  ~/reference/genome/human_g1k_v37/human_g1k_v37.fasta 1>human_g1k_v37.bwa_index.log 2>&1   &
mkdir -p ~/biosoft/GATK/resources/bundle/b37
cd ~/biosoft/GATK/resources/bundle/b37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.idx.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx.gz
gunzip 1000G_phase1.indels.b37.vcf.idx.gz
gunzip 1000G_phase1.indels.b37.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.b37.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.b37.vcf.idx.gz
mkdir -p ~/annotation/variation/human/dbSNP 
cd ~/annotation/variation/human/dbSNP 
## https://www.ncbi.nlm.nih.gov/projects/SNP/
## ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/
## ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/
nohup wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz &
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz.tbi

#数据分析--fastq2bam
module load java/1.8.0_91
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
INDEX=/home/jianmingzeng/reference/index/bwa/human_g1k_v37
GATK=/home/jianmingzeng/biosoft/GATK/GenomeAnalysisTK.jar
PICARD=/home/jianmingzeng/biosoft/picardtools/2.9.2/picard.jar
DBSNP=/home/jianmingzeng/annotation/variation/human/dbSNP/All_20160601.vcf.gz
SNP=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
INDEL=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
TMPDIR=/home/jianmingzeng/tmp/software
## samtools and bwa are in the environment 
## samtools Version: 1.3.1 (using htslib 1.3.1)
## bwa Version: 0.7.15-r1140
: '
'
## please keep the confige in three columns format, which are fq1 fq2 sampe
cat $1 |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    sample=${arr[2]}
    #####################################################
    ################ Step 1 : Alignment #################
    #####################################################
    echo bwa `date`
    bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 > $sample.sam
    echo bwa `date`
    #####################################################
    ################ Step 2: Sort and Index #############
    #####################################################
    echo SortSam `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$sample.sam OUTPUT=$sample.bam 
    samtools index $sample.bam
    echo SortSam `date`
    #####################################################
    ################ Step 3: Basic Statistics ###########
    #####################################################
    echo stats `date`
    samtools flagstat $sample.bam > ${sample}.alignment.flagstat
    samtools stats  $sample.bam > ${sample}.alignment.stat
    echo plot-bamstats -p ${sample}_QC  ${sample}.alignment.stat
    echo stats `date`
    #####################################################
    ####### Step 4: multiple filtering for bam files ####
    #####################################################
    ###MarkDuplicates###
    echo MarkDuplicates `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD MarkDuplicates \
    INPUT=$sample.bam OUTPUT=${sample}_marked.bam METRICS_FILE=$sample.metrics  
    echo MarkDuplicates `date`
    ###FixMateInfo###
    echo FixMateInfo `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD FixMateInformation \
    INPUT=${sample}_marked.bam OUTPUT=${sample}_marked_fixed.bam SO=coordinate  
    samtools index ${sample}_marked_fixed.bam
    echo FixMateInfo `date`
    echo ${sample}_marked_fixed.bam >>files.bamlist
    rm $sample.sam $sample.bam ${sample}_marked.bam
done 
samtools merge -@ 5  -b files.bamlist  merged.bam
samtools index merged.bam

#GATK重新处理bam文件
module load java/1.8.0_91
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
INDEX=/home/jianmingzeng/reference/index/bwa/human_g1k_v37
GATK=/home/jianmingzeng/biosoft/GATK/GenomeAnalysisTK.jar
PICARD=/home/jianmingzeng/biosoft/picardtools/2.9.2/picard.jar
DBSNP=/home/jianmingzeng/annotation/variation/human/dbSNP/All_20160601.vcf.gz
KG_SNP=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
Mills_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
KG_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.indels.b37.vcf
TMPDIR=/home/jianmingzeng/tmp/software
## samtools and bwa are in the environment 
## samtools Version: 1.3.1 (using htslib 1.3.1)
## bwa Version: 0.7.15-r1140
sample='merge'
###RealignerTargetCreator###
echo RealignerTargetCreator `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T RealignerTargetCreator \
-I ${sample}.bam -R $GENOME -o ${sample}_target.intervals \
-known $Mills_indels -known $KG_indels -nt 5
echo RealignerTargetCreator `date`
###IndelRealigner###
echo IndelRealigner `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T IndelRealigner \
-I ${sample}.bam -R $GENOME -targetIntervals ${sample}_target.intervals \
-o ${sample}_realigned.bam -known $Mills_indels -known $KG_indels 
echo IndelRealigner `date`
###BaseRecalibrator###
echo BaseRecalibrator `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T BaseRecalibrator \
-I ${sample}_realigned.bam -R $GENOME -o ${sample}_temp.table -knownSites $DBSNP
echo BaseRecalibrator `date`
###PrintReads###
echo PrintReads `date`
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T PrintReads \
-R $GENOME -I ${sample}_realigned.bam -o ${sample}_recal.bam -BQSR ${sample}_temp.table
samtools index ${sample}_recal.bam
echo PrintReads `date`
###delete_intermediate_files###

#variant calling by gatk hc
module load java/1.8.0_91
GENOME=/home/jianmingzeng/reference/genome/human_g1k_v37/human_g1k_v37.fasta
INDEX=/home/jianmingzeng/reference/index/bwa/human_g1k_v37
GATK=/home/jianmingzeng/biosoft/GATK/GenomeAnalysisTK.jar
PICARD=/home/jianmingzeng/biosoft/picardtools/2.9.2/picard.jar
DBSNP=/home/jianmingzeng/annotation/variation/human/dbSNP/All_20160601.vcf.gz
KG_SNP=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
Mills_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
KG_indels=/home/jianmingzeng/biosoft/GATK/resources/bundle/b37/1000G_phase1.indels.b37.vcf
TMPDIR=/home/jianmingzeng/tmp/software
## samtools and bwa are in the environment 
## samtools Version: 1.3.1 (using htslib 1.3.1)
## bwa Version: 0.7.15-r1140
fq1=P_jmzeng_DHG09057_AH33KVALXX_L1_1.clean.fq.gz
fq2=P_jmzeng_DHG09057_AH33KVALXX_L1_2.clean.fq.gz
sample='merge'
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T HaplotypeCaller  \
-R $GENOME -I ${sample}_recal.bam --dbsnp $DBSNP  \
-stand_emit_conf 10 -o  ${sample}_recal_raw.snps.indels.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xmx40g -jar $GATK -T HaplotypeCaller  \
-R $GENOME -I ${sample}_realigned.bam --dbsnp $DBSNP  \
-stand_emit_conf 10 -o  ${sample}_realigned_raw.snps.indels.vcf

