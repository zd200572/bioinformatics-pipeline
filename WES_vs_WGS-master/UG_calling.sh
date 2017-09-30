#!/bin/bash	

# Aziz Belkadi July 2014 
# WES and WGS SNV calling process using Unified Genotyper
# This script generates a WES SNV vcf file and a WGS SNV vcf file per sample
# Needs GATK Tools in the repertory (/resources)
# The WES bam files should be in (/WES/bam/)
# The WGS bam files should be in (/WGS/bam/)








mkdir /WES/vcf
mkdir /WGS/vcf


######################################
############ Declaration #############
######################################


reference=/resources/human_g1k_v37.fasta	
sureselect_intervals=/resources/sureselect/SureSelect-v4plusUTR_baits.interval_list	
gatk=/resources/apps/gatk/GenomeAnalysisTKLite-3.1-1/GenomeAnalysisTKLite.jar	
dbsnp=/resources/dbnsp/dbSNP135-All.vcf





##########################################################################
############ Calling a combined vcf file using GATK pipeline #############
##########################################################################


## List all WES bam files in the /WES repertory

ls /WES/bam/*.bam > /WES/bam/bam.list
ls /WGS/bam/*.bam > /WGS/bam/bam.list



## Calling WES in a combined vcf file

java -Xmx20g \
-jar $gatk \
-nt 8 \
--dbsnp $dbsnp \
-T UnifiedGenotyper \
-l INFO \
-rf BadCigar \
-L $intervals \
--interval_padding 50 \
-R $reference \
-I /WES/bam/bam.list \
-o /WES/vcf/combined.vcf



## Split WES combined vcf file across individuals

START=1
END=$(ls /WES/bam/*.bam | wc -l)-1

for (( c=$START; c<=$END; c++ ))
do
ind=$((9+$c))
l=`grep "CHROM" /WES/vcf/combined.vcf | cut -f$ind`
java -Xmx10g \
-jar $gatk \
-R $reference \
-T SelectVariants \
--variant /WES/vcf/combined.vcf \
-o /WES/vcf/$l.vcf \
-sn $l \
-env
done



## Calling WGS in a combined vcf file

java -Xmx20g \
-jar $gatk \
-nt 8 \
--dbsnp $dbsnp \
-T UnifiedGenotyper \
-l INFO \
-rf BadCigar \
-L $intervals \
--interval_padding 50 \
-R $reference \
-I /WGS/bam/bam.list \
-o /WGS/vcf/combined.vcf



## Split WGS combined vcf file across individuals

START=1
END=$(ls /WGS/bam/*.bam | wc -l)-1

for (( c=$START; c<=$END; c++ ))
do
ind=$((9+$c))
l=`grep "CHROM" /WGS/vcf/combined.vcf | cut -f$ind`
java -Xmx10g \
-jar $gatk \
-R $reference \
-T SelectVariants \
--variant /WGS/vcf/combined.vcf \
-o /WGS/vcf/$l.vcf \
-sn $l \
-env
done

