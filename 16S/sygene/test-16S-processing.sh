fq1=$1
fq2=$2
dirname = $3
mkdir -p $3 && cd $3
mkdir -p fastq-join & mkdir -p temp
mkdir -p result
#cd fastq-join
#拼接
join_paired_ends.py -f /media/vd1/171213_16s/$1 -r /media/vd1/171213_16s/$2 -o fastq-join