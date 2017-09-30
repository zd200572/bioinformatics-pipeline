## annovar and GATK 
 
 
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
 
mkdir -p  ~/biosoft/myBin
echo 'export PATH=/home/jianmingzeng/biosoft/myBin/bin:$PATH' >>~/.bashrc 
source ~/.bashrc
cd ~/biosoft
mkdir cmake &&  cd cmake
wget http://cmake.org/files/v3.3/cmake-3.3.2.tar.gz
tar xvfz cmake-3.3.2.tar.gz
cd cmake-3.3.2 
./configure --prefix=/home/jianmingzeng/biosoft/myBin
make
make install 
 
 
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
 
# ## Download and install tabix 
# cd ~/biosoft
# mkdir tabix &&  cd tabix
# # http://genometoolbox.blogspot.com/2013/11/installing-tabix-on-unix.html
# tar xvfj tabix-0.2.6.tar.bz2 
# cd tabix-0.2.6
# make
 
# cd ~/biosoft
# ## http://samtools.github.io/bcftools/
# mkdir htslib &&  cd htslib  
# git clone git://github.com/samtools/htslib.git 
# cd htslib
# make 
 
 
## Download and install bcftools
## http://www.htslib.org/download/
## http://www.htslib.org/doc/bcftools-1.0.html
cd ~/biosoft
mkdir bcftools &&  cd bcftools
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar xvfj bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1 
make
cp bcftools /home/jianmingzeng/biosoft/myBin
~/biosoft/myBin/bin/bcftools --help
 
 
## Download and install vcftools
## https://vcftools.github.io/index.html 
## http://vcftools.sourceforge.net/specs.html
cd ~/biosoft
mkdir vcftools &&  cd vcftools
 
# wget https://codeload.github.com/vcftools/vcftools/legacy.zip/master -O  vcftools-vcftools-v0.1.14-24-gac1bfd5.zip 
# unzip vcftools-vcftools-v0.1.14-24-gac1bfd5.zip 
# mv vcftools-vcftools-ac1bfd5 vcftools-v0.1.14-24
# cd vcftools-v0.1.14-24
# export PERL5LIB=/home/jianmingzeng/biosoft/vcftools/vcftools-v0.1.14-24/src/perl/
# ./autogen.sh 
# ./configure     --prefix=/home/jianmingzeng/biosoft/myBin
# make 
# make install 
# ~/biosoft/myBin/bin/vcftools --help
 
wget https://sourceforge.net/projects/vcftools/files/vcftools_0.1.13.tar.gz 
tar zxvf vcftools_0.1.13.tar.gz
cd  vcftools_0.1.13
make
 
 
 
## Download and install ANNOVAR  
cd ~/biosoft
# The latest version of ANNOVAR (2016Feb01) can be downloaded here (registration required)
# http://www.openbioinformatics.org/annovar/annovar_download_form.php 
mkdir ANNOVAR  &&  cd ANNOVAR  
 
  
 
 
 
## Download and install samstat
## http://samstat.sourceforge.net/
## http://www.htslib.org/doc/samtools.html
cd ~/biosoft
mkdir samstat &&  cd samstat
wget http://liquidtelecom.dl.sourceforge.net/project/samstat/samstat-1.5.tar.gz
tar zxvf  samstat-1.5.tar.gz 
cd samstat-1.5 
./configure --prefix=/home/jianmingzeng/biosoft/myBin
make
make install
~/biosoft/myBin/bin/samstat --help
 
 
## Download and install picardtools
## https://sourceforge.net/projects/picard/
## https://github.com/broadinstitute/picard
cd ~/biosoft
mkdir picardtools &&  cd picardtools
wget http://ncu.dl.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip
unzip picard-tools-1.119.zip 
 
 
## Download and install freebayes
## https://github.com/ekg/freebayes
## http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html
cd ~/biosoft
mkdir freebayes &&  cd freebayes
## wget -O freebayes-master.zip  https://codeload.github.com/ekg/freebayes/zip/master
## unzip freebayes-master.zip
wget http://clavius.bc.edu/~erik/freebayes/freebayes-5d5b8ac0.tar.gz
tar xzvf freebayes-5d5b8ac0.tar.gz
cd freebayes
make
 ~/biosoft/freebayes/freebayes/bin/freebayes
  
cd ~/biosoft
## https://sourceforge.net/projects/varscan/files/
## http://varscan.sourceforge.net/index.html
mkdir VarScan  &&  cd VarScan  
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar 
 
cd ~/biosoft
mkdir SnpEff &&  cd SnpEff
##  http://snpeff.sourceforge.net/
##  http://snpeff.sourceforge.net/SnpSift.html 
##  http://snpeff.sourceforge.net/SnpEff_manual.html
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip 
## java -jar snpEff.jar download GRCh37.75
## java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.75 example.vcf > example_snpeff.vcf
unzip snpEff_latest_core.zip
 
## https://github.com/najoshi/sickle
cd ~/biosoft
mkdir sickle && cd sickle
wget https://codeload.github.com/najoshi/sickle/zip/master -O sickle.zip
unzip sickle.zip
cd sickle-master
make
~/biosoft/sickle/sickle-master/sickle -h
 
cd ~/biosoft
## http://www.usadellab.org/cms/?page=trimmomatic
## http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
mkdir Trimmomatic && cd Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip 
unzip Trimmomatic-0.36.zip 
java -jar ~/biosoft/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar -h
 
 
 
## Download and install bedtools
cd ~/biosoft
mkdir bedtools &&  cd bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
## Length: 19581105 (19M) [application/octet-stream] 
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
#~/biosoft/bedtools/bedtools2/bin/
 
 
## Download and install PeakRanger
cd ~/biosoft
mkdir PeakRanger &&  cd PeakRanger
wget https://sourceforge.net/projects/ranger/files/PeakRanger-1.18-Linux-x86_64.zip 
## Length: 1517587 (1.4M) [application/octet-stream]
unzip PeakRanger-1.18-Linux-x86_64.zip
~/biosoft/PeakRanger/bin/peakranger -h
 
## Download and install bowtie
cd ~/biosoft
mkdir bowtie &&  cd bowtie
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip 
#Length: 27073243 (26M) [application/octet-stream]
#Saving to: "download"   ## I made a mistake here for downloading the bowtie2 
mv download  bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
 
## Download and install BWA
cd ~/biosoft
mkdir bwa &&  cd bwa
#http://sourceforge.net/projects/bio-bwa/files/
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2 
tar xvfj bwa-0.7.15.tar.bz2 # x extracts, v is verbose (details of what it is doing), f skips prompting for each individual file, and j tells it to unzip .bz2 files
cd bwa-0.7.15
make
#export PATH=$PATH:/path/to/bwa-0.7.15 # Add bwa to your PATH by editing ~/.bashrc file (or .bash_profile or .profile file)
# /path/to/ is an placeholder. Replace with real path to BWA on your machine
#source ~/.bashrc
   
 
 
 
## Download and install macs2  
## // https://pypi.python.org/pypi/MACS2/
cd ~/biosoft
mkdir macs2 &&  cd macs2
## just stick to PyPI release: https://pypi.python.org/pypi/MACS2
wget https://pypi.python.org/packages/9f/99/a8ac96b357f6b0a6f559fe0f5a81bcae12b98579551620ce07c5183aee2c/MACS2-2.1.1.20160309.tar.gz -O MACS2-2.1.1.tar.gz 
tar zxvf  MACS2-2.1.1.tar.gz 
cd  MACS2-2.1.1.20160309/
python setup.py install --user 
## https://docs.python.org/2/install/
~/.local/bin/macs2 --help
#wget https://codeload.github.com/taoliu/MACS/zip/master -O MACS-master.zip
#unzip MACS-master.zip
#cd MACS-master 
## So you can't just pull github snapshot then run setup.py to install MACS2. Instead
 
# ImageMagick
cd ~/biosoft
mkdir ImageMagick &&  cd ImageMagick
## http://www.imagemagick.org/ 
 
 
cd ~/biosoft
mkdir weblogo &&  cd weblogo
## http://weblogo.berkeley.edu/
wget http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz
tar zxvf weblogo.2.8.2.tar.gz
cd weblogo
export PATH=$PATH:/home/jianmingzeng/biosoft/weblogo/weblogo
source ~/.bashrc
 
cd ~/biosoft
mkdir Ghostscript &&  cd Ghostscript
# http://www.ghostscript.com/download/gsdnld.html
# http://www.ghostscript.com/doc/9.20/Readme.htm
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs920/ghostscript-9.20-linux-x86_64.tgz 
tar zxvf ghostscript-9.20-linux-x86_64.tgz
cp ghostscript-9.20-linux-x86_64/gs-920-linux_x86_64  ~/biosoft/myBin/bin/gs
## make sure the "gs" program is executable 
 
## Download and install homer (Hypergeometric Optimization of Motif EnRichment)
## // http://homer.salk.edu/homer/
## // http://blog.qiubio.com:8080/archives/3024 
## The commands gs, seqlogo, blat, and samtools should now work from the command line
cd ~/biosoft
mkdir homer &&  cd homer
wget http://homer.salk.edu/homer/configureHomer.pl 
perl configureHomer.pl -install
perl configureHomer.pl -install hg19
 
## Download and install SWEMBL
cd ~/biosoft
mkdir SWEMBL &&  cd SWEMBL
#### readme: http://www.ebi.ac.uk/~swilder/SWEMBL/beginners.html
wget http://www.ebi.ac.uk/~swilder/SWEMBL/SWEMBL.3.3.1.tar.bz2 
tar xvfj SWEMBL.3.3.1.tar.bz2
cd SWEMBL.3.3.1/
make
## error 
 
## Download and install SISSRs
cd ~/biosoft
mkdir SISSRs &&  cd SISSRs
#### readme: https://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/SISSRs-Manual.pdf
wget http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/sissrs_v1.4.tar.gz
tar xzvf sissrs_v1.4.tar.gz
~/biosoft/SISSRs/sissrs.pl
 
## Download and install SISSRs
cd ~/biosoft
mkdir QuEST &&  cd QuEST
#### http://mendel.stanford.edu/SidowLab/downloads/quest/
wget http://mendel.stanford.edu/SidowLab/downloads/quest/QuEST_2.4.tar.gz
tar xzvf QuEST_2.4.tar.gz
cd QuEST_2.4 
 
 
## Download and install fastqc
cd ~/biosoft
mkdir fastqc &&  cd fastqc
wget http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
 
 
## Download and install CEAS    
## http://liulab.dfci.harvard.edu/CEAS/download.html
 
cd ~/biosoft
mkdir CEAS  &&  cd CEAS 
wget  http://liulab.dfci.harvard.edu/CEAS/src/CEAS-Package-1.0.2.tar.gz
tar zxvf CEAS-Package-1.0.2.tar.gz
cd  CEAS-Package-1.0.2
python setup.py install --user 
## http://liulab.dfci.harvard.edu/CEAS/usermanual.html
 ~/.local/bin/ceas --help  
mkdir annotation  &&  cd annotation  
wget http://liulab.dfci.harvard.edu/CEAS/src/hg19.refGene.gz ; gunzip hg19.refGene.gz 
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz ##  gunzip refGene.txt.gz ; mv refGene.txt  hg19refGene.txt
 
## Download and install CEAS    
## http://liulab.dfci.harvard.edu/CEAS/download.html
 
cd ~/biosoft
mkdir crossmap  &&  cd crossmap 
pip install CrossMap --user
## http://crossmap.sourceforge.net/#use-pip-to-install-crossmap
 
mkdir chain_files  &&  cd chain_files  
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz 
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz ##  gunzip refGene.txt.gz ; mv refGene.txt  hg19refGene.txt
 
# Usage: CrossMap.py bed ~/biosoft/crossmap/chain_files/mm9ToMm10.over.chain.gz  test.mm9.bed3
 
 
 
cd ~/biosoft
# http://www.broadinstitute.org/cancer/cga/rnaseqc_run
# http://www.broadinstitute.org/cancer/cga/rnaseqc_download
mkdir RNA-SeQC  &&  cd RNA-SeQC 
#### readme: http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/RNA-SeQC_Help_v1.1.2.pdf
wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v1.1.8.jar 
 
 
#TopHat+Cufflinks+ pipeline
 
## Download and install TopHat 
# https://ccb.jhu.edu/software/tophat/index.shtml
cd ~/biosoft
mkdir TopHat  &&  cd TopHat 
#### readme: https://ccb.jhu.edu/software/tophat/manual.shtml
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar xzvf tophat-2.1.1.Linux_x86_64.tar.gz 
ln -s tophat-2.1.1.Linux_x86_64 current 
 
# ~/biosoft/TopHat/current/tophat2
 
## Download and install Cufflinks 
#  http://cole-trapnell-lab.github.io/cufflinks/
cd ~/biosoft
mkdir Cufflinks  &&  cd Cufflinks 
#### readme: http://cole-trapnell-lab.github.io/cufflinks/manual/
#### install:http://cole-trapnell-lab.github.io/cufflinks/install/
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xzvf cufflinks-2.2.1.Linux_x86_64.tar.gz 
ln -s cufflinks-2.2.1.Linux_x86_64 current
~/biosoft/Cufflinks/current/cufflinks
 
#HISAT-Stringtie2-Ballgown pipeline
 
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
 
## Download and install StringTie
## https://ccb.jhu.edu/software/stringtie/  ## https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
cd ~/biosoft
mkdir StringTie &&  cd StringTie 
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.3.Linux_x86_64.tar.gz 
tar zxvf  stringtie-1.2.3.Linux_x86_64.tar.gz
ln -s stringtie-1.2.3.Linux_x86_64 current
# ~/biosoft/StringTie/current/stringtie
 
cd ~/biosoft
mkdir RSEM &&  cd RSEM 
wget https://codeload.github.com/deweylab/RSEM/tar.gz/v1.2.31
mv v1.2.31  RSEM.v1.2.31.tar.gz 
tar zxvf RSEM.v1.2.31.tar.gz  
 
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
 
 
## Download and install kallisto
## https://pachterlab.github.io/kallisto/starting
cd ~/biosoft
mkdir kallisto &&  cd kallisto 
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz
#tar zxvf  
 
## Download and install Sailfish
## http://www.cs.cmu.edu/~ckingsf/software/sailfish/  ## 
cd ~/biosoft
mkdir Sailfish &&  cd Sailfish 
wget   https://github.com/kingsfordgroup/sailfish/releases/download/v0.9.2/SailfishBeta-0.9.2_DebianSqueeze.tar.gz 
#tar zxvf  
 
## Download and install salmon
## http://salmon.readthedocs.io/en/latest/salmon.html ## 
cd ~/biosoft
mkdir salmon &&  cd salmon 
## https://github.com/COMBINE-lab/salmon
#tar zxvf  
 
 
cd ~/biosoft
mkdir GDC  &&  cd GDC  
# https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
# http://gdc-docs.nci.nih.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/
wget https://gdc.nci.nih.gov/files/public/file/gdc-client_v1.2.0_Ubuntu14.04_x64.zip 
unzip gdc-client_v1.2.0_Ubuntu14.04_x64.zip
 
 
 
cd ~/biosoft/myBin/bin
## http://hgdownload.cse.ucsc.edu/admin/exe/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedSort
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer
 
## Download and install variationtoolkit
## https://code.google.com/archive/p/variationtoolkit/downloads 
cd ~/biosoft
mkdir variationtoolkit &&  cd variationtoolkit  
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/variationtoolkit/archive.tar.gz
tar zxvf archive.tar.gz 
cd variationtoolkit
make
 
## Download and install transvar
# http://bioinformatics.mdanderson.org/main/Transvar
cd ~/biosoft
# https://bitbucket.org/wanding/transvar
mkdir transvar &&  cd transvar  
wget https://bitbucket.org/wanding/transvar/get/v2.1.19.20160221.zip 
unzip v2.1.19.20160221.zip 
cd wanding-transvar-5dd8a7366999 
python setup.py install --user 
 
cd ~/biosoft
# http://kobas.cbi.pku.edu.cn/download.php 
mkdir kobas &&  cd kobas  
# wget http://kobas.cbi.pku.edu.cn/download_file.php?type=seq_pep&filename=hsa.pep.fasta.gz 
# wget http://kobas.cbi.pku.edu.cn/download_file.php?type=sqlite3&filename=hsa.db.gz 
wget http://kobas.cbi.pku.edu.cn/kobas-2.1.1/kobas-2.1.1.tar.gz 
tar zxvf kobas-2.1.1.tar.gz 
cd kobas-2.1.1_20160822
 
# * Download the KOBAS organism data package (organism.db.gz) from KOBAS Backend databases download website
# * Download the KOBAS specific species data package from KOBAS Backend databases download website (for example, hsa.db.gz)
# * Download the specific sequence file from KOBAS sequence files download website (for example, hsa.pep.fasta.gz)
# * `gunzip organism.db.gz`
# * `gunzip hsa.db.gz`
# * Move all databases into ${kobas_home}/sqlite3/ (for example, organism.db, hsa.db)
# * `gunzip hsa.pep.fasta.gz`
# * Move the fasta sequence file into ${kobas_home}/seq_pep/
# * `makeblastdb -in hsa.pep.fasta -dbtype prot`
 
         
pip install RPy2 --user 
pip install Numpy --user 
pip install Pandas --user 
pip install BioPython --user 
pip install matplotlib --user 
pip install PySQLite --user 
 
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
 
 
 
 
## Download and install bamtools
## https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
cd ~/biosoft
mkdir bamtools &&  cd bamtools  
git clone git://github.com/pezmaster31/bamtools.git 
cd bamtools
cmake --version  ## BamTools requires CMake (version >= 2.6.4).
mkdir build &&  cd build 
cmake ../ 
make
~/biosoft/bamtools/bamtools/bin/bamtools
 
## Download and install BAMStats
## http://bamstats.sourceforge.net/
## https://sourceforge.net/projects/bamstats/files/
cd ~/biosoft
mkdir BAMStats &&  cd BAMStats  
wget https://nchc.dl.sourceforge.net/project/bamstats/BAMStats-1.25.zip 
unzip BAMStats-1.25.zip
  
#java -jar  ~/biosoft/BAMStats/BAMStats-1.25/BAMStats-1.25.jar  --help
 
## Download and install Qualimap 
## http://qualimap.bioinfo.cipf.es/
cd ~/biosoft
mkdir Qualimap &&  cd Qualimap  
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip 
## readme  http://qualimap.bioinfo.cipf.es/doc_html/index.html
## example results :http://kokonech.github.io/qualimap/HG00096.chr20_bamqc/qualimapReport.html 
unzip qualimap_v2.2.1.zip 
~/biosoft/bamtools/bamtools/bin/bamtools
~/biosoft/Qualimap/qualimap_v2.2.1/qualimap --help ## --java-mem-size=4G
 
## modify ~/.bashrc by adding PATH=$PATH:~/.local/bin/
 
 
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
 
 
 
 
 
 
cd ~/biosoft
mkdir ngsplot &&  cd ngsplot  
## download by yourself :https://drive.google.com/drive/folders/0B1PVLadG_dCKN1liNFY0MVM1Ulk  
tar -zxvf ngsplot-2.61.tar.gz
tar zxvf ngsplot.eg.bam.tar.gz
 
echo 'export PATH=/home/jianmingzeng/biosoft/ngsplot/ngsplot/bin:$PATH' >>~/.bashrc  
echo 'export NGSPLOT=/home/jianmingzeng/biosoft/ngsplot/ngsplot' >>~/.bashrc 
source ~/.bashrc
 
install.packages("doMC", dep=T)
install.packages("caTools", dep=T)
install.packages("utils", dep=T)
  
source("http://bioconductor.org/biocLite.R")
biocLite( "BSgenome" )
biocLite( "Rsamtools" )
biocLite( "ShortRead" )
 
cd ~/biosoft
mkdir breakdancer &&  cd breakdancer  
# http://breakdancer.sourceforge.net/
# you need to install 2 perl module by yourself : http://breakdancer.sourceforge.net/moreperl.html
 
wget https://sourceforge.net/projects/breakdancer/files/breakdancer-1.1.2_2013_03_08.zip 
unzip breakdancer-1.1.2_2013_03_08.zip 
cd breakdancer-1.1.2/cpp
make  ##something wrong !
 
 
 
## usage: http://breakdancer.sourceforge.net/pipeline.html
 
cd ~/biosoft
# http://boevalab.com/FREEC/
mkdir Control-FREEC && cd Control-FREEC
# https://github.com/BoevaLab/FREEC/releases
wget https://github.com/BoevaLab/FREEC/archive/v10.3.zip 
unzip v10.3.zip 
# https://www.ncbi.nlm.nih.gov/pubmed/22155870
# http://boevalab.com/FREEC/tutorial.html
# http://samtools.sourceforge.net/pileup.shtml
 
 
cd ~/biosoft
# https://github.com/dellytools/delly
mkdir delly && cd delly 
# git clone --recursive https://github.com/dellytools/delly.git
# cd delly 
# make all
# make PARALLEL=1 -B src/delly
wget https://github.com/dellytools/delly/releases/download/v0.7.6/delly_v0.7.6_linux_x86_64bit 
chmod 777 delly_v0.7.6_linux_x86_64bit 
~/biosoft/delly/delly_v0.7.6_linux_x86_64bit  --help 
## delly call -t DEL -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam
## ./delly/src/bcftools/bcftools view delly.bcf > delly.vcf
## The SV type can be DEL, DUP, INV, TRA, or INS for deletions, tandem duplications, inversions, translocations and small insertions, respectively.
## In addition, you can filter input reads more stringently using -q 20 and -s 15.
 
 
 
 
cd ~/biosoft
# https://www.cog-genomics.org/plink2/data#merge_list
mkdir PLINK && cd PLINK 
wget https://www.cog-genomics.org/static/bin/plink170113/plink_linux_x86_64.zip 
unzip plink_linux_x86_64.zip
~/biosoft/PLINK/plink
 
## Download and install Scalpel
cd ~/biosoft
mkdir Scalpel &&  cd Scalpel
wget https://downloads.sourceforge.net/project/scalpel/scalpel-0.5.3.tar.gz  
tar zxvf scalpel-0.5.3.tar.gz
cd scalpel-0.5.3
make
~/biosoft/Scalpel/scalpel-0.5.3/scalpel-discovery  --help
~/biosoft/Scalpel/scalpel-0.5.3/scalpel-export  --help
 
cd ~/biosoft
# https://www.cog-genomics.org/plink2/data#merge_list
mkdir firehose && cd firehose 
wget http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip
unzip firehose_get_latest.zip 
~/biosoft/firehose/firehose_get
~/biosoft/firehose/firehose_get -tasks clinical analyses latest brca 
 
cd ~/biosoft
# https://www.cog-genomics.org/plink2/data#merge_list
mkdir fastpop && cd fastpop 
wget https://sourceforge.net/projects/fastpop/files/FastPop.tar.gz 
wget https://jaist.dl.sourceforge.net/project/fastpop/FastPop_Instruction.pdf
tar zxvf FastPop.tar.gz  
 
pip install cnvkit --user