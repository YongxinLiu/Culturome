#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
input=result/isolate_sample.txt
# 输出目录
output=result/purity
# 文库名称
l=L1

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    plot_purity.sh
Version:     1.0
Date:        2020/4/20
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: plot purity in heatmap 
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Jingying Zhang, Yong-Xin Liu, et. al. 
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).
-------------------------------------------------------------------------------
Version 1.0 2020/4/3

OPTIONS:
	-i directory of Culturome, default result/isolate_sample.txt
	-o output directory, default result/purity
	-d library ID, default L1
	-? show help of script 显示帮助

Example:
plot_purity.sh -i ${input} -o ${output} -l ${l}

EOF
}


# 参数解析 Analysis parameter
while getopts "i:o:l:" OPTION
do
	case $OPTION in
		i)
			input=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		l)
			l=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

# Each well purity in plate
# Fetch purity in result table
mkdir -p ${output}
awk 'NR==FNR{a[$1]=$3}NR>FNR{print $1"\t"a[$1]}' ${input} ${l}.txt | tail -n+2 | sed 's/\t$/\t0/' \
    > ${output}/${l}.txt
# Format list into plate format
format_list2plate.pl -i ${output}/${l}.txt \
    -o ${output}/${l}/

# Batch plot each plate
list=`ls ${output}/${l}/|cut -f 1 -d '.'|cut -f 2 -d 'P'|sort|uniq`
time for plate in $list;do 
    plot_pheatmap.sh -i ${output}/${l}/${l}P${plate}.plate \
        -o ${output}/${l}/${l}P${plate}.png
done

echo "Purify heatmap in ${output}/${l}/*.png"
