wget  ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.20110907.txt
cat CCDS.20110907.txt |perl -alne '{/[(.*?)]/;next unless $1;$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$" foreach split/,/,$exons;}' >hg19exon.bed

cat snp.vcf | java -jar  ~/biosoft/SnpEff/snpEff/SnpSift.jar intervals hg19exon.bed >hg19exon.snp.vcf
cat indel.vcf | java -jar  ~/biosoft/SnpEff/snpEff/SnpSift.jar intervals hg19exon.bed >hg19exon.indel.vcf

cat hg19_exon.snp.vcf |grep -v "^#" |cut -f 3 |grep '\.' |wc

cat hg19_exon.snp.vcf |perl -alne '{print if $F[2] eq "."}' |less -S
cat hg19_exon.indel.vcf |perl -alne '{print if $F[2] eq "."}' |less -S
