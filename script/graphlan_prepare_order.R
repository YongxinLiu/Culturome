#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <yxliu@genetics.ac.cn>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：统计ASV表为样本列表和非冗余ASV列表
# Functions: Identify non-redundancy isolates


options(warn = -1) # Turn off warning

# 1.2 参数 Parameters #----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="ASV table in reads count [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy_8.txt",
                help="metadata file or metadata [default %default]"),
    make_option(c("-a", "--abundance"), type="numeric", default="0.1",
                help="Mean relative abundance in percentage [default %default]"),
    make_option(c("-n", "--number"), type="numeric", default="100",
                help="Top N taxonomy to plot [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/graphlan",
                help="Output quantile value for filter feature table [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print(opts)
dir.create(opts$output)

# Version 1.0, Based on ASV table and taxonomy, output graphlan
# 版本 1.0, 基于ASV表和7级物种注释文件，输出物种树图

# 参数设置
# 按丰度筛选，如0.01即代表0.01%，即万分之一
# abundance = 0.01
# 按数量筛选，如150即代表最高丰度的150个特征
# number = 150

# 读取输入文件 #----
otutab = read.table(opts$input, sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F, comment.char = "")
taxonomy = read.table(opts$taxonomy, sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F, comment.char = "")


# 数据筛选#----
# 标准化并求均值
norm = as.data.frame(t(t(otutab)/colSums(otutab,na=T)*100))
# 丰度由大到小排序
idx = order(rowMeans(norm), decreasing = T)
norm = norm[idx,]
# 按丰度筛选
idx = rowMeans(norm) > opts$abundance
filtered_otutab = norm[idx,]
# 按数量筛选
filtered_otutab = head(filtered_otutab, opts$number)
# 添加均值并保留4位小数
filtered_otutab = round(cbind(rowMeans(filtered_otutab), filtered_otutab), digits = 4)
colnames(filtered_otutab)[1] = "Mean"
filtered_otutab$Mean = filtered_otutab$Mean / max(filtered_otutab$Mean)
# 对应过滤物种注释
idx = rownames(filtered_otutab) %in% rownames(taxonomy)
filtered_otutab = filtered_otutab[idx,]
filtered_taxonomy = taxonomy[rownames(filtered_otutab),]

# 保存输出文件
# 过滤的OTU表
write.table("OTUID\t", file=paste0(opts$output, "/filtered_otutab.txt"), append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
suppressWarnings(write.table(filtered_otutab, file=paste0(opts$output, "/filtered_otutab.txt"), append = T, sep="\t", quote=F, row.names=T, col.names=T))
# 过滤的物种注释
write.table("OTUID\t", file=paste0(opts$output, "/filtered_taxonomy.txt"), append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
suppressWarnings(write.table(filtered_taxonomy, file=paste0(opts$output, "/filtered_taxonomy.txt"), append = T, sep="\t", quote=F, row.names=T, col.names=T))

### 绘制树骨架 #----

# 输入文件为筛选后的taxonomy文件：filtered_taxonomy.txt

# 读取筛选后的文件，不设置行名
tax = read.table(paste0(opts$output, "/filtered_taxonomy.txt"), sep="\t", header = TRUE, stringsAsFactors = F)
# 筛选门-属5级+OTUID
tree = data.frame(tax[,c(3:7,1)], stringsAsFactors = F)
# head(tree)
## clarify taxonomy，解决不同级别重名问题，为可识别级别，且与Greengene格式保持一致
tree[,1] = paste("p__",tree[,1],sep = "")
tree[,2] = paste("c__",tree[,2],sep = "")
# tree[,3] = paste("o__",tree[,3],sep = "")
tree[,4] = paste("f__",tree[,4],sep = "")
tree[,5] = paste("g__",tree[,5],sep = "")
# save tree backbone, 按点分隔格式

# 解决科标签重名问题
idx = tree[,3] %in% "Unassigned"
# 方法1. 重名标签添加数字编号，但结果有太多Unassigned
# tree[idx,4] = paste0(tree[idx,4], 1:length(tree[idx,4]))
# 方法2. 过滤掉科末注释的条目，数量会减少，但图片更美观
tree = tree[!idx,]
# 简化一些代_的不规则科名
tree[,3] = gsub('_\\w*',"",tree[,3])
write.table (tree, file=paste0(opts$output, "/tree1_backbone.txt"), sep=".", col.names=F, row.names=F, quote=F)

# 列出现在有门、纲、目、科、属，用于设置与门对应的背景色
Phylum = unique(tree[,1]) 
Class = unique(tree[,2])
Order = unique(tree[,3])
Family = unique(tree[,4])
Genus = unique(tree[,5])

# 筛选四大菌门中的科并按门着色
# 修改为目，则将tree的4列改为3列，Family改为Order
pro = tree[tree[,1]=="p__Proteobacteria",3]
act = tree[tree[,1]=="p__Actinobacteria",3] 
bac = tree[tree[,1]=="p__Bacteroidetes",3]
fir = tree[tree[,1]=="p__Firmicutes",3]

# 对每个科进行标签、文字旋转、按门注释背景色
# 也可调整为其它级别，如Order, Class或Genus
label_color = data.frame(stringsAsFactors = F)
for (element in Order)
{
  # element
  anno = data.frame(stringsAsFactors = F)
  anno[1,1] = element
  anno[1,2] = "annotation"
  anno[1,3] = "*"
  # 设置文字旋转90度
  anno[2,1] = element
  anno[2,2] = "annotation_rotation"
  anno[2,3] = "90"
  # 设置背景色，四大门各指定一种色，其它为灰色
  anno[3,1] = element
  anno[3,2] = "annotation_background_color" 
  
  if (element %in% pro)
  {
      anno[3,3] = "#85F29B"
  } else if (element %in% act)
  {
      anno[3,3] = "#F58D8D"   
  } else if (element %in% fir)
  {
      anno[3,3] = "#F7C875"  
  } else if (element %in% bac)
  {
      anno[3,3] = "#91DBF6"   
  } else {
      anno[3,3] = "grey"   
  }
  label_color = rbind(label_color,anno)
}
write.table(label_color, paste0(opts$output, "/tree2_label_color.txt"), sep = "\t", quote = F,col.names = F,row.names = F, na="")