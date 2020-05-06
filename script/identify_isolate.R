#!/usr/bin/env Rscript

# Copyright 2016-2020 Yuan Qin <yqin@genetics.ac.cn> Yong-Xin Liu <yxliu@genetics.ac.cn>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：统计ASV表为样本列表和非冗余ASV列表
# Functions: Identify non-redundancy isolates


options(warn = -1) # Turn off warning

# 1.2 参数 Parameters #----
# 设置清华源加速下载
# (Optional) Set up Tsinghua Mirror to speed up download
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析包是否安装，没安装则安装，然后加载
# Determine whether the command line parsing package is installed, install it if it is not installed, then load
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/ASV_table.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy_8.txt",
                help="metadata file or metadata [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/isolate",
                help="Output quantile value for filter feature table [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)

# Version 1.0, Based on ASV table and taxonomy, output well info (purity, counts and taxonomy) and candidate wells of non-redundancy ASV
# 版本 1.0, 基于ASV表和7级物种注释文件，输出每个孔的信息，筛选每个孔的信息()，以及非冗余ASV的修行孔，纯度优先，数据量其次排序的Top 5


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","grid","scales","dplyr")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("scales")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("grid")))

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

# 2.1 读文件 Load #----

# Load ASV table
# Public file 1. "result/ASV_table.txt"  raw reads count of each ASV in each sample
otu_table = read.delim(opts$input, row.names= 1,  header=T, sep="\t")

# Load ASV metadata
# Public file 2. "result/taxonomy_8.txt"  taxonomy for each ASV, tab seperated
taxonomy = read.delim("result/taxonomy_8.txt", row.names= 1,header=T, sep="\t")
taxonomy$Full=paste(taxonomy$Phylum,taxonomy$Class,taxonomy$Order,taxonomy$Family,taxonomy$Genus,taxonomy$Species,sep = ";")

# Extract only those samples in common between the two tables
idx = rownames(otu_table) %in% rownames(taxonomy)
otu_table = otu_table[idx,]
taxonomy = taxonomy[rownames(otu_table),]


# 2.2 稀疏曲线 Rarefraction curve #----

sample_rare  =  function(df, count_cutoff = 3, length = 30, rep = 30){  
  result = data.frame(sample = c(0), richness = c(0))
  count=df
  # Set otu table to binary for easy calculate richness
  count[count < count_cutoff] = 0
  count[count >= count_cutoff] = 1
  # Sample number
  n = dim(count)[2]
  # 设置X轴各取样数量
  x = unique(as.integer(seq(1, n, length.out = length)))
  for (i in x){
    m = choose(n, i)
    if (m > rep){m = rep}
    # loop list and calculate richness
    for(j in 1:m) {
      idx = sample(n, i)
      temp = count[, idx, drop = F]
      # subset non-zero row
      temp1=temp[rowSums(temp)>0, , drop=F]
      # row number is accumulative OTUs
      result = rbind(result, c(i, dim(temp1)[1]))
    }
  }
  # remove start 0,0
  result = result[-1,]
  # factor draw as each box
  result$sample = factor(result$sample, levels = unique(result$sample))
  return(result)
}
rare_box = sample_rare(otu_table, count_cutoff = 1, length = 30, rep = 30)
p = ggplot(rare_box,aes(x=sample,y=richness, fill=sample))+
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.5, size = 0.1)+main_theme+ xlab("Well numbers")+ylab("ASV number") + 
  theme(legend.position = "NA",axis.text.x = element_text(angle = 45,vjust=1, hjust=1))
p
ggsave(paste(opts$output, "_rare_curve.pdf", sep=""), p, width = 89*2, height = 59, unit = "mm")
ggsave(paste(opts$output, "_rare_curve.png", sep=""), p, width = 89*2, height = 59, unit = "mm")

print(paste("Rarefaction curve in ",opts$output, "_rare_curve.pdf/png", sep=""))

# 2.3 样品丰度和物种 Sample taxonomy and abundance #----

# 样品中ASV种类与丰度统计
# 统计每个孔中的前三个OTU
sample=dim(otu_table)[2]
sample_stat=as.data.frame(cbind(1:sample,colnames(otu_table)))
rownames(sample_stat)=sample_stat$V2

for(i in 1:sample) {
  temp=otu_table[order(otu_table[,i],decreasing = TRUE),c(i,1)]
  sample_stat[i,3]=temp[1,1]
  sample_stat[i,4]=round(temp[1,1]/sum(temp[,1])*100, digits = 2)
  sample_stat[i,5]=rownames(temp[1,])
  sample_stat[i,6]=taxonomy[rownames(temp[1,]),]$Full
  sample_stat[i,7]=temp[2,1]
  if (temp[2,1] > 0){
    sample_stat[i,8]=round(temp[2,1]/sum(temp[,1])*100, digits = 2)
    sample_stat[i,9]=rownames(temp[2,])
    sample_stat[i,10]=taxonomy[rownames(temp[2,]),]$Full
  }else{
    sample_stat[i,8]="NA"
    sample_stat[i,9]="NA"
    sample_stat[i,10]="NA"
  }
}
# 选菌原则：同一OTU，纯度优先、数据量第二，最多选前三；输出OTU, purity-count-
sample_stat=sample_stat[,-1]
colnames(sample_stat)=c("SampleID","Count","Purity","ASV","Taxonomy","Count2","Purity2","ASV2","Taxonomy2")
sample_stat$ID=as.numeric(gsub("ASV_", "", sample_stat$ASV))
sample_stat = arrange(sample_stat, ID, desc(Purity), desc(Count))
sample_stat = sample_stat[,-10]
suppressWarnings(write.table(sample_stat, file=paste(opts$output, "_well.txt", sep=""), append = T, sep="\t", quote=F, row.names=F, col.names=T))

print(paste("Well list in ",opts$output, "_well.txt", sep=""))

# 2.4 非冗余候选 Identify non-redundancy isolates #----

# 读取孔信息并选菌
sample_stat = read.delim(paste(opts$output, "_well.txt", sep=""), row.names= NULL, header=T, sep="\t", stringsAsFactors = F)

# 查看Top 1的ASV数量
print("Top1 ASV number:")
length(unique(sample_stat$ASV))
top1=unique(sample_stat$ASV)

# 将竖表整理为横表选菌，选Top 5，ID+纯度+丰度

## 目标菌建立列表，至少两列矩阵为数据框
otu_stat=as.data.frame(cbind(1:2,1:2))
sample=dim(sample_stat)[1]
otu="ASV_0"
# 每个ASV的候选孔数量初始值，最多5
j=0 
# ASV编号，从0起，一般200-500
m=0 
for(i in 1:sample) {
  # i=1
  if (otu != as.character( sample_stat[i,"ASV"])){
    otu=as.character( sample_stat[i,"ASV"])
    j=0
    m=m+1
    otu_stat[m,1]=as.character(sample_stat[i,"ASV"])
    otu_stat[m,2]=as.character(sample_stat[i,"Taxonomy"])
  }
  if (j==0){
    otu_stat[m,3]=as.character(sample_stat[i,"SampleID"])
    otu_stat[m,4]=as.numeric(sample_stat[i,"Purity"])
    otu_stat[m,5]=as.integer(sample_stat[i,"Count"])
  }
  if (j==1){
    otu_stat[m,6]=as.character(sample_stat[i,"SampleID"])
    otu_stat[m,7]=as.numeric(sample_stat[i,"Purity"])
    otu_stat[m,8]=as.integer(sample_stat[i,"Count"])
  }  
  if (j==2){
    otu_stat[m,9]=as.character(sample_stat[i,"SampleID"])
    otu_stat[m,10]=as.numeric(sample_stat[i,"Purity"])
    otu_stat[m,11]=as.integer(sample_stat[i,"Count"])
  }  
  if (j==3){
    otu_stat[m,12]=as.character(sample_stat[i,"SampleID"])
    otu_stat[m,13]=as.numeric(sample_stat[i,"Purity"])
    otu_stat[m,14]=as.integer(sample_stat[i,"Count"])
  }  
  if (j==4){
    otu_stat[m,15]=as.character(sample_stat[i,"SampleID"])
    otu_stat[m,16]=as.numeric(sample_stat[i,"Purity"])
    otu_stat[m,17]=as.integer(sample_stat[i,"Count"])
  }  
  j=j+1
}

## 删除ASV中非数字，并转换为数值
# otu_stat$V1=as.integer(gsub("ASV_","",otu_stat$V1,perl=TRUE))
## 按ASV数字顺序排序
# otu_stat = arrange(otu_stat, V1)
## 重命名行名
colnames(otu_stat)=c("ID","Taxonomy","Well1","Purity1","Count1","Well2","Purity2","Count2","Well3","Purity3","Count3","Well4","Purity4","Count4","Well5","Purity5","Count5")
## 写入文件
suppressWarnings(write.table(otu_stat, file=paste(opts$output, "_ASV.txt", sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=T))

# print(paste("Output in ", opts$output, "_ASV.txt", "finished. 选菌候选表", sep = ""))
# print(paste("Output in ", opts$output, "_well.txt", "finished. 每个孔的信息，包括丰度最高两种菌", sep = ""))

print(paste("ASV list in ",opts$output, "_ASV.txt", sep=""))


# 2.4 孔序列种类分布 Unique ASV distribution #----

# 统计每个样品的ASV数量
TopN = 7
# 转换为2元表
otu_table2 = otu_table
otu_table2[otu_table > 0] = 1
Number = as.data.frame(colSums(otu_table2))
# 设置图例数量
df = as.data.frame(table(Number), stringsAsFactors = F)
other = sum(df[TopN:dim(df)[1], ]$Freq)
dfN = head(df, n=TopN-1)
df = rbind(dfN, c(paste0(">=",TopN), other))
df$Number = factor(df$Number, levels = df$Number)
df$Freq = as.numeric(df$Freq)
# p = ggplot(df, aes(Number, Freq, fill = Number)) + 
#   geom_bar(stat="identity",width=0.7) 
# p
# ggsave(paste(opts$output, ".ASVdistribution.pdf", sep=""), p, width = 4, height = 3)
# ggsave(paste(opts$output, ".ASVdistribution.png", sep=""), p, width = 4, height = 3)

# 统计每个样品属的数量
genus_table = read.delim("result/tax/sum_g.txt", row.names= 1,  header=T, sep="\t")
genus_table2 = genus_table
genus_table2[genus_table > 0] = 1
Number = as.data.frame(colSums(genus_table2))
# 设置图例数量
df2 = as.data.frame(table(Number), stringsAsFactors = F)
other = sum(df2[TopN:dim(df2)[1], ]$Freq)
df2N = head(df2, n=TopN-1)
df2 = rbind(df2N, c(paste0(">=",TopN), other))
df2$Number = factor(df2$Number, levels = df2$Number)
df2$Freq = as.numeric(df2$Freq)
# p = ggplot(df2, aes(Number, Freq, fill = Number)) + 
#   geom_bar(stat="identity",width=0.7) 
# p
# ggsave(paste(opts$output, ".GenusDistribution.pdf", sep=""), p, width = 4, height = 3)
# ggsave(paste(opts$output, ".GenusDistribution.png", sep=""), p, width = 4, height = 3)

# 属和ASV混合
df$tax = "ASV"
df2$tax = "Genus"
dfAll = rbind(df,df2)
p = ggplot(dfAll, aes(Number, Freq, fill = tax)) + 
  geom_bar(stat="identity", position = "dodge",width=0.7) + main_theme +
  theme(legend.position = c(0.8, 0.8))
p
ggsave(paste(opts$output, ".Distribution.pdf", sep=""), p, width = 89, height = 59*0.7, unit = "mm")
ggsave(paste(opts$output, ".Distribution.png", sep=""), p, width = 89, height = 59, unit = "mm")

print(paste("Distribution of ASV/Genus number of all wells in ",opts$output, ".Distribution.pdf/png", sep=""))
