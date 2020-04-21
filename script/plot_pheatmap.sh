#!/bin/bash
set -e

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    plot_pheatmap.sh
Revision:    1.0
Date:        2017/10/18
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: Plot matrix heatmap by pheatmap
Notes:       Input matrix must have row and column name
-------------------------------------------------------------------------------
Version 1.0 2017/10/18

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, L1_split.count

# 1. input.txt, a matrix must have row and column name
Plate   1       2       3       4       5       6       7       8       9       10      11      12
A       881     66      35      13      14      3       1403    5       1402    104     0       0
B       745     1       0       660     523     0       1123    1155    0       866     8       41
C       1142    645     1412    32      1348    1193    31      779     1233    7       1094    557
D       14      1525    1083    34      1423    1155    1285    1138    15      4       0       598
E       6       0       1061    1052    1175    25      1196    11      1       0       552     1092
F       384     4       1306    0       1401    1       1021    789     0       963     561     327
G       998     620     24      1139    1       1100    9       861     0       0       11      4
H       0       2       1       0       605     112     0       8       0       4       0       0

# Output file
# 1. input.txt.pdf

OPTIONS:
	-e execuate Rscript, TRUE or FALSE, default TRUE
	-h figure height, default 5 inch
	-i input filename, default input.txt
	-o output filename, default input.txt.pdf
	-s text size, default 7
	-w figure width, default 8
	-I install package TRUE or FALSE, default false
	-c cell width, default 30
	-C cell height, default 30
	-S scale, row, col, none, default none


Example:
	plot_pheatmap.sh -i input.txt -o input.txt.pdf
EOF
}


# Default parameter
execute='TRUE'
height=2.7
input='input.txt'
output='input.pdf'
text_size=7
width=3.5
ist='FALSE' # install package, default FALSE
cellwidth=16
cellheight=16
scale='none'
cluster_rows=F
cluster_cols=F
labels_col='c(1:12)'
display_numbers=T
number_format='%.f'


# Analysis parameter
while getopts "e:h:i:o:s:w:I:c:C:S:" OPTION
do
	case $OPTION in
		e)
			execute=$OPTARG
			;;
		h)
			height=$OPTARG
			;;
		i)
			input=$OPTARG
			output="${input}.pdf"
			;;
		o)
			output=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		I)
			ist=$OPTARG
			;;
		c)
			cellwidth=$OPTARG
			;;
		C)
			cellheight=$OPTARG
			;;
		l)
			labels_col=$OPTARG
			;;
		L)
			labels_row=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done


# 输出R脚本
cat <<END >plot_pheatmap.r
# Install related packages
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2","pheatmap"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
#setwd("${output}")
library("ggplot2")
library("reshape2")
library("pheatmap")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=${text_size}),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=${text_size}),
                    text=element_text(family="sans", size=${text_size}))

# 绘制热图
data = read.table("${input}", header=T, row.names= 1,sep="\t") # dataframe
#pheatmap(data, cellwidth=${cellwidth}, cellheight=${cellwidth}, scale="${scale}", cluster_rows=${cluster_rows}, cluster_cols=${cluster_cols}, labels_col=${labels_col}, labels_row=labels_row, display_numbers=${display_numbers}, number_format="${number_format}", filename = "${output}", width=${width}, height=${height})
pheatmap(data, cellwidth=${cellwidth}, cellheight=${cellwidth}, scale="${scale}",
  cluster_rows=${cluster_rows}, cluster_cols=${cluster_cols}, 
  display_numbers=${display_numbers}, number_format="${number_format}",
  filename = "${output}", width=${width}, height=${height}, fontsize = ${text_size})

END


if test "${execute}" == "TRUE";
then
	Rscript plot_pheatmap.r
	rm plot_pheatmap.r
	rm -rf Rplots.pdf
fi
