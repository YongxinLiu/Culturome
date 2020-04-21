#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
input=temp/ASV.fa.tax
# 输出目录
output=result/
# 临时目录
database=result/ASV_table.txt


# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    tax_sum.sh
Version:     1.0
Date:        2020/4/6
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: Summary taxonomy in each level
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Jingying Zhang, Yong-Xin Liu, et. al. 
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).
-------------------------------------------------------------------------------
Version 1.0 2020/4/6

OPTIONS:
	-i directory of Culturome, default temp/ASV.fa.tax
	-o output directory, default result/taxonomy
	-d ASV table, default result/ASV_table.txt
	-? show help of script 显示帮助

Example:
tax_sum.sh -i ${input} -o ${output} -d ${database}

EOF
}


# 参数解析 Analysis parameter
while getopts "i:o:d:" OPTION
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
		?)
			usage
			exit 1
			;;
	esac
done

# Fill blank to unassianged
sed -i 's/\t$/\td:Unassigned/' ${input}

# Summary each taxonomic level
mkdir -p ${output}/tax
for i in p g;do \
	usearch -sintax_summary ${input} \
	-otutabin ${database} \
	-output ${output}/tax/sum_${i}.txt \
	-rank ${i} ; done
	
# Format adjustment
sed -i 's/(//g;s/)//g;s/\"//g;s/\/Chloroplast//g;s/\-/_/g;s/\//_/' ${output}/tax/sum_*.txt
cut -f 1,4 ${input} | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > ${output}/taxonomy_2.txt
awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
    ${output}/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
    > ${output}/taxonomy_8.txt
sed -i 's/#//g;s/ //g' ${output}/taxonomy_8.txt
echo "Result in ${output}taxonomy_8.txt"