#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
input=script/profiling_ASV.fa 
# 临时目录
database=script/profiling_ASVtab.txt
# 参考序列
reference=result/ASV.fa
# 输出目录
output=result/cross_reference.txt

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    cross_reference.sh
Version:     1.0
Date:        2020/5/6
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: Blast 16S profile sequences with cultured ASV, and calculated cultured rate
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Jingying Zhang, Yong-Xin Liu, et. al. 
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).
-------------------------------------------------------------------------------
Version 1.0 2020/5/6

OPTIONS:
	-i 16S representitive sequences, default script/profiling_ASV.fa
	-d 16S ASV table, default profiling_ASVtab.txt
	-r cultured ASV as reference, default result/ASV.fa
	-o ASV with cultured marker and relative abundance, default result/cross_reference.txt
	-? show help of script

Example:

cross_reference.sh \
    -i ${input} \
    -d ${database} \
    -r ${reference} \
    -o ${output}
EOF
}


# 参数解析 Analysis parameter
while getopts "i:o:d:r:" OPTION
do
	case $OPTION in
		i)
			input=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		d)
			database=$OPTARG
			;;
		r)
			reference=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done


# 比对
mkdir -p temp
makeblastdb -in ${reference} -dbtype nucl
grep  -c '>' ${input}
blastn -query ${input} \
    -db ${reference} \
    -out temp/profile.blastn \
    -outfmt 6 -num_alignments 1
wc -l temp/profile.blastn

# 筛选>97%的可培养ASV
awk '$3>=97' temp/profile.blastn > temp/profile.blastn97
wc -l temp/profile.blastn97

# 添加丰度信息
usearch -otutab_counts2freqs ${database} -output temp/profiling_ASVtab.txt
awk '{a=a+$2}END{print a}' temp/profiling_ASVtab.txt

# 统计可培养丰度
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]="Cultured"} NR>FNR {print $0,a[$1]}' temp/profile.blastn97 temp/profiling_ASVtab.txt > ${output}
echo "Cultured ASVs/OTUs number:"
grep 'Cultured' ${output} | wc -l
echo "Cultured relative abundance of ASVs/OTUs:"
grep 'Cultured' ${output} |  awk '{a=a+$2}END{print a}'
