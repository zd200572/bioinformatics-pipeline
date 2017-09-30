## all data come from : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81916
cut -f 3 config.txt |while read id ; do wget $id 2>/dev/null  ;done 

ls *sra |while read id; do  ~/biosoft/sratoolkit/sratoolkit.2.6.3-centos_linux64/bin/fastq-dump  --gzip --split-3  $id;done

ls *.fastq.gz |xargs ~/biosoft/fastqc/FastQC/fastqc -t 10

nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz  &

reference=/home/jianmingzeng/reference/index/hisat/mm10/genome 
~/biosoft/HISAT/current/hisat2  -p 5 -x $reference -1 SRR3589959_1.fastq.gz  -2  SRR3589959_2.fastq.gz -S control_1.sam  2>control_1.log
~/biosoft/HISAT/current/hisat2  -p 5 -x $reference -1 SRR3589960_1.fastq.gz  -2  SRR3589960_2.fastq.gz -S Akap95_1.sam  2>Akap95_1.log
~/biosoft/HISAT/current/hisat2  -p 5 -x $reference -1 SRR3589961_1.fastq.gz  -2  SRR3589961_2.fastq.gz -S control_2.sam  2>control_2.log
~/biosoft/HISAT/current/hisat2  -p 5 -x $reference -1 SRR3589962_1.fastq.gz  -2  SRR3589962_2.fastq.gz -S Akap95_2.sam  2>Akap95_2.log


ls *sam |while read id;do (nohup samtools sort -n -@ 5 -o ${id%%.*}.Nsort.bam $id   & );done

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz 
ls *.Nsort.bam |while read id;do (nohup samtools view  $id | ~/.local/bin/htseq-count -f sam -s no -i gene_name - ~/reference/gtf/gencode/gencode.vM12.annotation.gtf   1>${id%%.*}.geneCounts 2>${id%%.*}.HTseq.log   & );done

ls *sam |while read id;do (nohup samtools sort -@ 5 -o ${id%%.*}.Psort.bam $id    &);done

ls *.Psort.bam |while read id;do (nohup  samtools index  $id   & );done
ls *.Psort.bam |while read id;do (nohup bamCoverage -b $id -o  ${id%%.*}.bw    & );done

# you need create the config.txt file by yourself, below is the content: 
# GSM2177726	control_rep1	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX180/SRX1801303/SRR3589959/SRR3589959.sra
# GSM2177727	Akap95_rep1	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX180/SRX1801304/SRR3589960/SRR3589960.sra
# GSM2177728	control_rep2	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX180/SRX1801305/SRR3589961/SRR3589961.sra
# GSM2177729	Akap95_rep2	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX180/SRX1801306/SRR3589962/SRR3589962.sra


# you need to install some software by yourself,they are sratoolkit,fastqc,HISAT,samtools,deeptools,htseq-count 


## Download and install sratoolkit
## http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
## http://www.ncbi.nlm.nih.gov/books/NBK158900/
cd ~/biosoft
mkdir sratoolkit &&  cd sratoolkit
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.3/sratoolkit.2.6.3-centos_linux64.tar.gz
##
##  Length: 63453761 (61M) [application/x-gzip]
##  Saving to: "sratoolkit.2.6.3-centos_linux64.tar.gz"
tar zxvf sratoolkit.2.6.3-centos_linux64.tar.gz
~/biosoft/sratoolkit/sratoolkit.2.6.3-centos_linux64/bin/fastdump -h

## Download and install fastqc
cd ~/biosoft
mkdir fastqc &&  cd fastqc
wget http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip

## Download and install HISAT 
# https://ccb.jhu.edu/software/hisat2/index.shtml
cd ~/biosoft
mkdir HISAT  &&  cd HISAT 
#### readme: https://ccb.jhu.edu/software/hisat2/manual.shtml
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
ln -s hisat2-2.0.4  current 
## ~/biosoft/HISAT/current/hisat2-build
## ~/biosoft/HISAT/current/hisat2  

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


## Download and install deepTools
## https://github.com/fidelram/deepTools
## http://deeptools.readthedocs.io/en/latest/content/example_usage.html
pip install pyBigWig --user 

cd ~/biosoft
mkdir deepTools &&  cd deepTools  
git clone https://github.com/fidelram/deepTools ## 130M,
cd deepTools
python setup.py install --user
## 17 tools in ~/.local/bin/
~/.local/bin/deeptools 

## Download and install HTSeq  
## http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
## https://pypi.python.org/pypi/HTSeq
cd ~/biosoft
mkdir HTSeq &&  cd HTSeq
wget  https://pypi.python.org/packages/72/0f/566afae6c149762af301a19686cd5fd1876deb2b48d09546dbd5caebbb78/HTSeq-0.6.1.tar.gz 
tar zxvf HTSeq-0.6.1.tar.gz
cd HTSeq-0.6.1
python setup.py install --user 
~/.local/bin/htseq-count  --help
## ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/
## http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/liftOver/
## GRCm38/mm10 (Dec, 2011) 
## ls *bam |while read id;do ( ~/.local/bin/htseq-count  -f bam  $id   genecode/mm9/gencode.vM1.annotation.gtf.gz  1>${id%%.*}.gene.counts ) ;done 
## ls *bam |while read id;do ( ~/.local/bin/htseq-count  -f bam -i exon_id  $id   genecode/mm9/gencode.vM1.annotation.gtf.gz  1>${id%%.*}.exon.counts ) ;done

