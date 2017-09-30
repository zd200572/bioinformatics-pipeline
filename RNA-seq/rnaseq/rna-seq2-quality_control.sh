#rna-seq2-quality_control.sh
# quality control

cp /home/training/data/RNA-Seq/example2-* ~/proj2/reads/

cd ~/proj2/fastqc/

fastqc -f fastq -o ./ ../reads/example2-*

   

# back to workpath



cd ~/