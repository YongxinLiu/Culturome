#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
wd=~/github/Culturome
# 输出目录
output=result/graphlan/
# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    graphlan_plot.sh
Version:     1.0
Date:        2020/4/3
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: Based on usearch alph_div, draw alpha boxplot and statistics, also can draw any matrix to boxplot by selected group and feature
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
	-i directory of Culturome, default ~/github/Culturome
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
			wd=$OPTARG
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

# 配置文件
cfg=${wd}/script/cfg/


## 基本树绘图#----

rm -rf ${output}/track*
# 生成树的默认参数，可手动调整更多样式
cat ${cfg}/global.cfg ${output}/tree2_label_color.txt > ${output}/track0
# 合并所有的注释，接下来会生成更多track，使树更复杂
cat ${output}/track* > ${output}/graphlan_annotate.txt
## 注释树
#graphlan_annotate.py --annot ${output}/graphlan_annotate.txt \
#  ${output}/tree1_backbone.txt ${output}/graphlan.xml
## 绘图，size决定图片大小，越大字越小
#graphlan.py ${output}/graphlan.xml ${output}/graphlan1_tree.pdf --size 5


# 我们需要从树文件中获得节点名称，并添加注释数据。

#如获得结点的丰度，在下面很多注释都会基于丰度信息

# 获得最终出图的结点ID
cut -f 6 -d '.' ${output}/tree1_backbone.txt > ${output}/tree1_backbone.id
# 注释结果丰度均值
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' \
  ${output}/filtered_otutab.txt ${output}/tree1_backbone.id > ${output}/tree1_backbone.mean

## 形状标签有无

# 样式1. 如筛选丰度，用紫色方块标出大于千分之5的结点

# # 环1筛选千分之五的结果注释为方块，cfg/ring1.cfg中的m代表紫色，R代表方块
# cat ${cfg}/ring1.cfg <(awk '$2>0.5' ${output}/tree1_backbone.mean | \
#   cut -f 1 | sed 's/$/\tring_shape\t1\tR/') > ${output}/track1
# 
# # 绘图，加第一环矩形，展示丰度大于千万的特征
# cat ${output}/track* > ${output}/graphlan_annotate.txt
# graphlan_annotate.py --annot ${output}/graphlan_annotate.txt \
#   ${output}/tree1_backbone.txt ${output}/graphlan.xml
# graphlan.py ${output}/graphlan.xml ${output}/graphlan1_rectangle.pdf --size 5


# 样式2. 如筛选丰度，用第二环位置橙色倒三角标出小于千分之5的结点

# 注释：ring2.cfg为第二环，颜色y为yellow橙色，注释track中也为2

# # 环1筛选千分之五的结果注释为方块，cfg/ring1.cfg中的m代表紫色，R代表方块
# cat ${cfg}/ring2.cfg <(awk '$2<=0.5' ${output}/tree1_backbone.mean | \
#   cut -f 1 | sed 's/$/\tring_shape\t2\tv/') > ${output}/track2
# 
# # 绘图，加第一环矩形，展示丰度大于千万的特征
# cat ${output}/track* > ${output}/graphlan_annotate.txt
# graphlan_annotate.py --annot ${output}/graphlan_annotate.txt \
#   ${output}/tree1_backbone.txt ${output}/graphlan.xml
# graphlan.py ${output}/graphlan.xml ${output}/graphlan2_triangle.pdf --size 5


## 热图展示丰度

# 添加所有样品均值作为热图，作为第3环。本质上热图即环形条带的透明度

# 环1筛选千分之五的结果注释为方块，cfg/ring1.cfg中的m代表紫色，R代表方块
cat ${cfg}/heat3.cfg <(sed 's/\t/\tring_alpha\t3\t/g' ${output}/tree1_backbone.mean) > ${output}/track3

# 绘图，加第一环矩形，展示丰度大于千万的特征
cat ${output}/track* > ${output}/graphlan_annotate.txt
graphlan_annotate.py --annot ${output}/graphlan_annotate.txt ${output}/tree1_backbone.txt ${output}/graphlan.xml
graphlan.py ${output}/graphlan.xml ${output}/graphlan.pdf --size 5
# graphlan.py ${output}/graphlan.xml ${output}/graphlan.png --size 5


## 柱状图显示丰度

# 环1筛选千分之五的结果注释为方块，cfg/ring1.cfg中的m代表紫色，R代表方块
# cat ${cfg}/bar4.cfg <(sed 's/\t/\tring_height\t4\t/g' ${output}/tree1_backbone.mean) > ${output}/track4
# 
# # 绘图，加第一环矩形，展示丰度大于千万的特征
# cat ${output}/track* > ${output}/graphlan_annotate.txt
# graphlan_annotate.py --annot ${output}/graphlan_annotate.txt \
#   ${output}/tree1_backbone.txt ${output}/graphlan.xml
# graphlan.py ${output}/graphlan.xml ${output}/graphlan.pdf --size 5
# graphlan.py ${output}/graphlan.xml ${output}/graphlan.png --size 5
