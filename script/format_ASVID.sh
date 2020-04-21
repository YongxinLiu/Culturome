#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
input=temp/Zotus.fa
# 输出目录
output=result/otu.fa
# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    fomrat_ASVID.sh
Version:     1.0
Date:        2020/4/3
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: Based on usearch alph_div, draw alpha boxplot and statistics, also can draw any matrix to boxplot by selected group and feature
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
https://doi.org/10.1007/s11427-018-9284-4
-------------------------------------------------------------------------------
Version 1.0 2020/4/3

OPTIONS:
	-i directory of Culturome, default temp/Zotus.fa
	-o output directory, default result/graphlan
	-? show help of script 显示帮助

Example:
graphlan_plot.sh -i ${wd} -o ${output}

EOF
}


# 参数解析 Analysis parameter
while getopts "i:o:" OPTION
do
	case $OPTION in
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

awk 'BEGIN {n=1}; />/ {print ">ASV_" n; n++} !/>/ {print}' ${input} \
    > temp/otus.fa
cp temp/otus.fa ${output}
