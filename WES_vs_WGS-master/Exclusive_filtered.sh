#!/bin/bash	

# Aziz Belkadi July 2014 
# WES and WGS filtered exclusive vcf files per individual
# This script generates an exclusive WES SNV vcf file and an exclusive WGS SNV vcf file per individual
# Needs GATK Tools in the repertory (/resources)
# The WES filtered vcf files should be in (/Filter/WES)
# The WGS filtered vcf files should be in (/Filter/WGS)
# WES and WGS files for the same sample should have the same name








mkdir /Exclusive/Filter/WES
mkdir /Exclusive/Filter/WGS


######################################
############ Declaration #############
######################################


reference=/resources/human_g1k_v37.fasta	
gatk=/resources/apps/gatk/GenomeAnalysisTKLite-3.1-1/GenomeAnalysisTKLite.jar	




###########################################################
############ Select exclusive Variants in WES #############
###########################################################





## List all WES vcf files in the /Filter/WES repertory and all WES vcf files in the /Filter/WGS repertory

ls /Filter/WES/*.vcf > Filter/WES/vcf.list
ls /Filter/WGS/*.vcf > Filter/WGS/vcf.list



## Find discordant WES variants comparing to WGS

START=1
END=$(ls /Filter/WES/*.vcf | wc -l)

for (( c=$START; c<=$END; c++ ))
do
ind_WES=`head -$c /Filter/WES/vcf.list | tail -n 1`
ind_WGS=`head -$c /Filter/WGS/vcf.list | tail -n 1`
java -Xmx10g \
-jar $gatk \
-R $reference \
-T SelectVariants \
--variant "$ind_WES" \
--discordance "$ind_WGS" \
-o Exclusive/"$ind_WES"
done





## Find discordant WGS variants comparing to WES

START=1
END=$(ls /Filter/WGS/*.vcf | wc -l)

for (( c=$START; c<=$END; c++ ))
do
ind_WES=`head -$c /Filter/WES/vcf.list | tail -n 1`
ind_WGS=`head -$c /Filter/WGS/vcf.list | tail -n 1`
java -Xmx10g \
-jar $gatk \
-R $reference \
-T SelectVariants \
--variant "$ind_WGS" \
--discordance "$ind_WES" \
-o Exclusive/"$ind_WGS"
done

