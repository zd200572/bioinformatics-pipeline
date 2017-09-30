# generate folder

cd ~/

mkdir proj2

cd proj2

mkdir reads fastqc ref tophat

   

# copy tools

cd ~/

cp -fr /home/training/tools/bowtie ./tools/

cp -fr /home/training/tools/tophat ./tools/

cp -fr /home/training/tools/cufflinks ./tools/

   

# change path

export PATH=$PATH:~/bin/

   

# install bowtie

cd ~/tools/bowtie/

unzip bowtie2-2.0.5-linux-x86_64.zip

cd bowtie2-2.0.5/

cp bowtie2* ~/bin/

   

# install tophat

cd ~/tools/tophat/

tar -zxvf tophat-2.0.8.Linux_x86_64.tar.gz

cd tophat-2.0.8.Linux_x86_64/

cp * ~/bin/

   

# install cufflinks

cd ~/tools/cufflinks/

tar -zxvf cufflinks-2.0.2.Linux_x86_64.tar.gz

cd cufflinks-2.0.2.Linux_x86_64/

cp * ~/bin/

   

# back to workpath

echo "preparation and install complete..."

cd ~/