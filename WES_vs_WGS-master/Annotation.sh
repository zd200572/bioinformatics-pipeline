#!/bin/bash	

# Aziz Belkadi july 2014 
# Annotation
# Annovar software uziped in /resources/annovar repertory
# result in /Annotation/Exclusive/Filter/
# input exclusive and filtered vcf files in /Exclusive/Filter/



	







mkdir -p /Annotation/Exclusive/Filter/WES
mkdir -p /Annotation/Exclusive/Filter/WGS



######################################
############ Declaration #############
######################################

annovar_dir=/resources/annovar





#####################################
############ Annotation #############
#####################################



## List all WES vcf files in the /Exclusive/Filter/WES repertory and all WES vcf files in the /Exclusive/Filter/WGS repertory

ls /Exclusive/Filter/WES/*.vcf > /Exclusive/Filter/WES/vcf.list
ls /Exclusive/Filter/WGS/*.vcf > /Exclusive/Filter/WGS/vcf.list


## Annotation of WES variants

START=1
END=$(ls /Exclusive/Filter/WES/*.vcf | wc -l)

for (( c=$START; c<=$END; c++ ))
do
ind_WES=`head -$c /Exclusive/Filter/WES/vcf.list | tail -n 1`

## CONVERT VCF TO ANNOVAR FORMAT
$annovar_dir/convert2annovar.pl $ind_WES -format vcf4 > /Annotation/$ind_WES.variantCalls.vcf.annovar

## SUMMARIZE ANNOVAR
$annovar_dir/summarize_annovar.pl --buildver hg19 --verdbsnp 135 --ver1000g 1000g2012apr --veresp 6500 --remove /Annotation/$ind_WES.variantCalls.vcf.annovar $annovar_dir/humandb/

done






## Annotation of WGS variants

START=1
END=$(ls /Exclusive/Filter/WGS/*.vcf | wc -l)

for (( c=$START; c<=$END; c++ ))
do
ind_WGS=`head -$c /Exclusive/Filter/WGS/vcf.list | tail -n 1`

## CONVERT VCF TO ANNOVAR FORMAT
$annovar_dir/convert2annovar.pl $ind_WGS -format vcf4 > /Annotation/$ind_WGS.variantCalls.vcf.annovar

## SUMMARIZE ANNOVAR
$annovar_dir/summarize_annovar.pl --buildver hg19 --verdbsnp 135 --ver1000g 1000g2012apr --veresp 6500 --remove /Annotation/$ind_WGS.variantCalls.vcf.annovar $annovar_dir/humandb/

done


