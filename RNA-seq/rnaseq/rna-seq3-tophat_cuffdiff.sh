#rna-seq3-tophat_cuffdiff.sh
# index genome

cp /home/training/data/RNA-Seq/ref2.fa ~/proj2/ref/

cp /home/training/data/RNA-Seq/ann2.gtf ~/proj2/ref/

   

cd ~/proj2/ref/

bowtie2-build ref2.fa ref2

   

# generate bam file

cd ~/proj2/tophat/

tophat2 -o E2-1-thout ../ref/ref2 ../reads/example2-1.L.fq ../reads/example2-1.R.fq

tophat2 -o E2-2-thout ../ref/ref2 ../reads/example2-2.L.fq ../reads/example2-2.R.fq

   

# generate gtf file

cd ~/proj2/tophat/

cufflinks -o E2-1-clout E2-1-thout/accepted_hits.bam

cufflinks -o E2-2-clout E2-2-thout/accepted_hits.bam

   

# generate assemblies.txt

touch assemblies.txt

echo "./E2-1-clout/transcripts.gtf" >> assemblies.txt

echo "./E2-2-clout/transcripts.gtf" >> assemblies.txt

cd ~/proj2/tophat/

cuffmerge -s ../ref/ref2.fa assemblies.txt

   

# find differential genes based on merged.gtf

cuffdiff -o diff_out1 -b ../ref/ref2.fa -L E2-1,E2-2 -u merged_asm/merged.gtf ./E2-1-thout/accepted_hits.bam ./E2-2-thout/accepted_hits.bam

   

# find differential genes based on ann2.gtf

cuffdiff -o diff_out2 -b ../ref/ref2.fa -L E2-1,E2-2 -u ../ref/ann2.gtf ./E2-1-thout/accepted_hits.bam ./E2-2-thout/accepted_hits.bam