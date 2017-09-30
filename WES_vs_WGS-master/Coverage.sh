#!/bin/bash	

# Aziz Belkadi July 2014 
# WES and WGS compute coverage across exons
# Needs GATK Tools in the repertory (/resources)
# The WES bam files should be in (/WES/bam/)
# The WGS bam files should be in (/WGS/bam/)
# The refGene file in /resources repertory








mkdir /WES/coverage
mkdir /WGS/coverage


######################################
############ Declaration #############
######################################


reference=/resources/human_g1k_v37.fasta	
sureselect_intervals=/resources/sureselect/SureSelect-v4plusUTR_baits.interval_list	
gatk=/resources/apps/gatk/GenomeAnalysisTKLite-3.1-1/GenomeAnalysisTKLite.jar	
refGene=/resources/refgene_b37.txt






##########################################################################
############ Calling a combined vcf file using GATK pipeline #############
##########################################################################


## List all WES bam files in the /WES repertory

ls /WES/bam/*.bam > /WES/bam/bam.list
ls /WGS/bam/*.bam > /WGS/bam/bam.list



## WES coverage

java -Xmx10g \
-jar $gatk \
-nt 8 \
-T DepthOfCoverage \
-L $intervals \	
--interval_padding 50 \
-R $reference \
-I /WES/bam/bam.list \
-o /WES/coverage/WES_coverage \
-geneList $refGene \
-ct 4 \
-omitBaseOutput




## WES coverage

java -Xmx10g \
-jar $gatk \
-nt 8 \
-T DepthOfCoverage \
-L $intervals \	
--interval_padding 50 \
-R $reference \
-I /WGS/bam/bam.list \
-o /WGS/coverage/WGS_coverage \
-geneList $refGene \
-ct 4 \
-omitBaseOutput
