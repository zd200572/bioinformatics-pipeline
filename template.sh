#!/bin/bash
set -e

# 设置程序参数的缺省值，少用参数即可运行
# Default parameter
input=input.txt
output=output.txt
database=database.txt
execute='TRUE'

# 程序功能描述，每次必改程序名、版本、时间等；版本更新要记录清楚，一般可用-h/-?来显示这部分
# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    template.sh
Revision:    1.0
Date:        2017/6/24
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is solve parameter read and default
Notes:       Function of this script
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/6/24

# 输入输出文件格式和示例，非常有用，不知道格式怎么准备文件呀
# Input files: input.txt, can inclue many file

# 1. input.txt, design of expriment
SampleID    BarcodeSequence    group    gene    batch    description
WT.1    TAGCTT    WT    ggps9.10    2    double mutant of ggps9-ggps10, cause A/B down
WT.2    GGCTAC    WT    ggps9.10    2    double mutant of ggps9-ggps10, cause A/B down
WT.3    CGCGCG    WT    ggps9.10    2    double mutant of ggps9-ggps10, cause A/B down

# 2. database.txt, annotation of gene
ID    description
AT3G48300    Transcript factor

# Output file
1. Annotated samples & DE genes
Samples    ID    description
Wt    AT3G48300    Transcript factor

2. Volcano plot: vol_otu_SampleAvsSampleB.pdf

# 参数描述，写清功能的缺省值
OPTIONS:
    -d database file, default database.txt
    -i input file, recommend must give
    -o output file or output directory, default output.txt
    -h/? show help of script
Example:
    template.sh -i input.txt -d database.txt -o result.txt
EOF
}

# 解释命令行参数，是不是很面熟，其实是调用了perl语言的getopts包，
# Analysis parameter
while getopts "d:h:i:o:" OPTION
do
    case $OPTION in
        d)
            database=$OPTARG
            ;;
        h)
            usage
            exit 1
            ;;
        i)
            input=$OPTARG
            ;;
        o)
            output=$OPTARG
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done