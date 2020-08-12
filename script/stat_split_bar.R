#!/usr/bin/env Rscript

# Copyright 2016-2020 Yuan Qin <yqin@genetics.ac.cn> Yong-Xin Liu <yxliu@genetics.ac.cn>

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
    make_option(c("-i", "--input"), type="character", default="result/metadata.txt",
                help="Metadat [default %default]"),
    make_option(c("-d", "--database"), type="character", default="result/split/L1.txt",
                help="Reads counts of library [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/split",
                help="Output directory [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

print("You are using the following parameters:")
print(opts)

# Version 1.0, Plotting reads count distribution in barplot 
# 版本 1.0, 绘制每个孔的分布，和整合高到低的分布

library("ggplot2")
library("stringr")

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

# Public file 1. "design.txt"  Design of experiment
design = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char='') 

# Public file 2. "L1.txt"  raw reads count of each OTU in each sample
sample_count = read.table(paste(opts$database, sep=""), header=F, sep="\t") # dataframe
colnames(sample_count)=c("Sample","Count")

design$group = design$plate
design = design[,c("group", "plate", "library","well")]
df = cbind(sample_count, design[match(sample_count$Sample, rownames(design)), ]) 

# 提取孔的行号，用于着色
split1 = as.data.frame(str_split_fixed(df$well, "",2))
# split1 = split1[,6:7]
colnames(split1)[1] = "row"
colnames(split1)[2] = "num"

# split2 = as.data.frame(str_split_fixed(df$Sample, "",6))
# split2 = split2[,5:6]
# colnames(split2)[1] = "x"
# colnames(split2)[2] = "well"
# 
# split = cbind(split2,split1)
# split <- split[,2:3]
df_split <- cbind(df, split1)

df_split$well  = factor(df_split$well, levels=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12"))

# 2.1 绘制每板的分布#----
# 共48板分4行，每行12板，每板中按8行着色
p = ggplot(df_split, aes(x=well, y = Count, fill=row))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  xlab("Library")+ylab("Pair Reads count")+
  facet_wrap(~ group, scales = "fixed",  nrow=4)+
  theme(axis.text.x = element_text(angle = 90))+ 
  scale_fill_brewer(palette="Set1")+
  main_theme+labs(title = opts$database)
# p
ggsave(paste(opts$output, "/", basename(opts$database), ".well.pdf", sep=""), p, width = 16, height = 9)
ggsave(paste(opts$output, "/", basename(opts$database), ".well.png", sep=""), p, width = 16, height = 9)
print(paste("Well counts in ",opts$output, basename(opts$database), ".well.pdf/png", sep=""))

# 2.2 绘制整体的分布#----
p2 = ggplot(df_split, aes(x=Count, fill = "group"))+ 
  geom_histogram(binwidth = 50)+ 
  xlab("Count of reads")+ylab("Numbers")+
  theme(axis.text.x = element_text(angle = 90))+
  main_theme+labs(title=opts$database)+theme(legend.position="none")
# p2
ggsave(paste(opts$output, "/", basename(opts$database), ".histogram.pdf", sep=""), p2, width = 8, height = 5)
ggsave(paste(opts$output, "/", basename(opts$database), ".histogram.png", sep=""), p2, width = 8, height = 5)

print(paste("Histogram of well counts in ",opts$output, basename(opts$database), ".histogram.pdf/png", sep=""))

# 2.3 绘制每板的数据量#----
suppressWarnings(suppressMessages(library(dplyr)))
df = df[,c("Count","group")]
df$group = gsub("\\w+P","",df$group, perl = T)
plate_count = as.data.frame(df %>% group_by(group) %>% summarise_all(sum))
# 去除开始的0
# 转换为数值坐标变为1，10间隔刻度了
# plate_count$group = as.numeric(gsub("^0","",plate_count$group))
plate_count$group = gsub("^0","",plate_count$group)
plate_count$group = factor(plate_count$group, levels = plate_count$group)

p3 = ggplot(plate_count, aes(x=group, y=Count, fill=group))+ 
  geom_bar(stat="identity",position="dodge", width=0.7)+ 
  xlab("Plates")+ylab("Count of reads")+
  main_theme+theme(axis.text.x = element_text(angle = 90,vjust=1, hjust=1), legend.position="none") # + labs(title = opts$database)
# 查看绘图的数据和颜色 
# ggplot_build(p3)$data 
# 修改颜色，循环彩虹色
# p3 + scale_fill_manual(values = rep(rainbow(7),7)[1:length(plate_count$group)])
# 添加均值线
p3 = p3 + geom_hline(yintercept=mean(plate_count$Count), linetype=2, color="blue") 
# p3  
ggsave(paste(opts$output, "/", basename(opts$database), ".plate.pdf", sep=""), p3, width = 89, height = 59, unit = "mm")
ggsave(paste(opts$output, "/", basename(opts$database), ".plate.png", sep=""), p3, width = 8, height = 5)

write.table(plate_count, file=paste(opts$output, "/", basename(opts$database), ".plate.txt", sep=""),append = F, quote = F, sep="\t", eol = "\n", row.names = F, col.names = T)

print(paste("Plate counts in ",opts$output, basename(opts$database), ".plate.pdf/png", sep=""))

print("Analysis completed!!!")
