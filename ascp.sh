vol1/fastq/ERR262/ERR262997/ERR262997_1.fastq.gz

anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/sra/SRX/SRX248/SRX248556/SRR771547/SRR771547.sra

ascp -QT -l 6M -k1 -i asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR262/ERR262997/ERR262997_1.fastq.gz .

ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR262/ERR262997

ascp -QT -l 6M -k1 -i asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/sra/SRX/SRX248/SRX248556/SRR771547/SRR771547.sra

 ascp -i asperaweb_id_dsa.openssh -Tr -Q -l 6M -P33001 -L- -k1 era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR262/ERR262997/ERR262997_1.fastq.gz H:\
 
 ascp -QT -l 6M -k1 -i asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP017/SRP017311/SRR620204/SRR620204.sra d:\
 ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP017/SRP017311/SRR620204/SRR620204.sra
 
 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam
 
ascp -i asperaweb_id_dsa.openssh -Tr -Q -l 6M -P33001 -L- -k1 fasp-g1k@fasp.1000genomes.ebi.ac.uk:/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam