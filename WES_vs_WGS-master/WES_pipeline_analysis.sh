#!/bin/bash	

# Aziz Belkadi July 2014 
# Paired ends WES sequencing pipeline
# This script generates a final alignment bam file
# Maximum Exact Matching (MEM) alignment algorithm
# Needs BWA, GATK, Picard Tools in the same repertory (/resources)
# The two fastq files should have the same name with R1.fastq extention for the first end read file (Ex: Ind1_R1.fastq.gz) and R2.fastq for the second end read file (Ex: Ind1_R2.fastq.gz)
# The two fastq files should be in the same repertory (/WES/fastq)
# Take one parameter : Fastq files name without exetention 	


	
!/bin/bash

if [ $# -ne 1 ]
then
	echo "Usage: <SampleID>"
	exit 65
fi



############################
master_time_start=`date +%s`
############################

mkdir -p /WES/bam/$1-b37/logs
mkdir /WES/bam/$1-b37/tmpdir




######################################
############ Declaration #############
######################################

data_dir=/WES/fastq
out_dir=/WES/bam/$1-b37

reference=/resources/human_g1k_v37.fasta	
sureselect_intervals=/resources/sureselect/SureSelect-v4plusUTR_baits.interval_list	
bwa_dir=/resources/bwa-0.7.8	
picard_dir=/resources/apps/picard-tools-1.113	
gatk=/resources/apps/gatk/GenomeAnalysisTKLite-3.1-1/GenomeAnalysisTKLite.jar	
kg_mills=/resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.gz
kg_indels=/resources/1000G_phase1.indels.b37.vcf.gz
dbsnp=/resources/dbnsp/dbSNP135-All.vcf




########################################################
############ Alignment using GATK pipeline #############
########################################################


	
#BWA ALIGNMENT	
	
$bwa_dir/bwa mem \	
-t 8 \
-M \
$reference \
$data_dir/$1*R1*fastq*gz \
$data_dir/$1*R2*fastq*gz \
> $out_dir/$1.aln.sampe.sam


## SORT BAM AND ADD INFO
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx20g \
-jar $picard_dir/AddOrReplaceReadGroups.jar \
I=$out_dir/$1.aln.sampe.sam \
O=$out_dir/$1.aln.sampe.sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
RGID=$1 \
RGLB="pe" \
RGPU="HiSeq-2000" \
RGSM=$1 \
RGCN="Human Genetics of Infectious Disease" \
RGDS=$intervals--GRCh37 \
RGPL=illumina \
VALIDATION_STRINGENCY=SILENT \
&>> $out_dir/logs/$1.log

## CREATE INTERVALS FOR LOCAL REALIGNMENT
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx20g \
-jar $gatk \
-R $reference \
-L $intervals \
--interval_padding 50 \
-T RealignerTargetCreator \
-rf BadCigar \
-known $kg_mills \
-known $kg_indels \
-nt 1 \
-I $out_dir/$1.aln.sampe.sorted.bam \
-o $out_dir/$1.aln.sampe.sorted.bam.forRealigner.intervals \
--allow_potentially_misencoded_quality_scores \
&>> $out_dir/logs/$1.log

## PERFORM LOCAL REALIGNMENT
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx20g \
-jar $gatk \
-R $reference \
-T IndelRealigner \
-rf BadCigar \
-known $kg_mills \
-known $kg_indels \
--consensusDeterminationModel USE_READS \
-compress 0 \
-targetIntervals $out_dir/$1.aln.sampe.sorted.bam.forRealigner.intervals \
-I $out_dir/$1.aln.sampe.sorted.bam \
-o $out_dir/$1.aln.sampe.sorted.realigned.bam \
--allow_potentially_misencoded_quality_scores \
&>> $out_dir/logs/$1.log

## MARK DUPLICATES
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx20g \
-jar $picard_dir/MarkDuplicates.jar \
REMOVE_DUPLICATES=false \
M=$out_dir/metrics/$1.duplicate.metrics \
I=$out_dir/$1.aln.sampe.sorted.realigned.bam \
O=$out_dir/$1.aln.sampe.sorted.realigned.dedup.bam \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
&>> $out_dir/logs/$1.log

## BASE RECALIBRATOR
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx20g \
-jar $gatk \
-L $intervals \
--interval_padding 50 \
-R $reference \
--disable_indel_quals \
-T BaseRecalibrator \
-rf BadCigar \
-knownSites $dbsnp \
-knownSites $kg_mills \
-knownSites $kg_indels \
-I $out_dir/$1.aln.sampe.sorted.realigned.dedup.bam \
-o $out_dir/$1.recal_data.grp \
&>> $out_dir/logs/$1.log

## APPLY RECALIBRATION (PRINT READS)
java -Djava.io.tmpdir=$out_dir/tmpdir \
-Xmx10g \
-jar $gatk \
-R $reference \
-T PrintReads \
-rf BadCigar \
-BQSR $out_dir/$1.recal_data.grp \
-I $out_dir/$1.aln.sampe.sorted.realigned.dedup.bam \
-o $out_dir/$1.aln.sampe.sorted.realigned.dedup.recal.bam \
--allow_potentially_misencoded_quality_scores \
&>> $out_dir/logs/$1.log

	
	
	
#####################################
###### Remove temporary files #######
#####################################

rm $out_dir/$1.recal_data.grp
rm $out_dir/$1.aln.sampe.sorted.realigned.dedup.bai
rm $out_dir/$1.aln.sampe.sorted.realigned.dedup.bam
rm $out_dir/$1.aln.sampe.sorted.realigned.bai
rm $out_dir/$1.aln.sampe.sorted.realigned.bam
rm $out_dir/$1.aln.sampe.sorted.bam.forRealigner.intervals
rm $out_dir/$1.aln.sampe.sorted.bai
rm $out_dir/$1.aln.sampe.sorted.bam
rm $out_dir/$1.aln.sampe.sam
rm -r $out_dir/tmpdir

############################
master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))`; echo "$1 analysis completed in $master_time_exec seconds") >> $out_dir/logs/$1.log
############################
