#!/usr/bin/env Rscript

# Copyright 2016-2020 Yuan Qin <yqin@genetics.ac.cn> Yong-Xin Liu <yxliu@genetics.ac.cn>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. 
# NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
# Nature Biotechnology 37, 676-684, 
# doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# Clean enviroment object
rm(list=ls()) 

# 1.1 程序功能描述和主要步骤 #----

# 程序功能：筛选OTU表前，根据每板阴性对照reads值设置阈值
# Functions: OTU filter threshold

options(warn = -1) # Turn off warning

# 1.2 解析命令行 #----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="temp/otutab.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="result/metadata.txt",
                help="metadata file or metadata [default %default]"),
    make_option(c("-n", "--negative"), type="character", default="A12",
                help="Well ID of negative col [default %default]"),
    make_option(c("-t", "--threshold"), type="numeric", default="1",
                help="Quantile of reads [default %default]"),
    make_option(c("-p", "--positive"), type="character", default="B12",
                help="Posotive control well ID [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/fdr.txt",
                help="Output quantile value for filter feature table [default %default]") 
)
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}

print("You are using the following parameters:")
print(opts)

library(ggplot2)
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

# 2. 依赖关系检查、安装和加载#----

# suppressWarnings(suppressMessages(library(amplicon)))

# 3. 读取输入文件#----

# Public file 1. "metadata.txt"  metadata of experiment
metadata = read.table(opts$metadata, header=T, row.names= 1, sep="\t", comment.char='') 
metadata$SampleID = rownames(metadata)
metadata = metadata[,c("SampleID","library","plate","well")]

# Public file 2. "otutab.txt"
otutab = read.delim(opts$input, row.names= 1,  header=T, sep="\t")

# 交叉筛选OTU表和实验设计共有的样品
otutab_f = otutab[,colnames(otutab) %in% rownames(metadata)]
metadata_f = metadata[rownames(metadata) %in% colnames(otutab),]

readsload= as.data.frame(colSums(otutab_f))
colnames(readsload) = "readssum"

# metadata$group = metadata$plate

## 3.1 补充降噪过程中去掉的孔号#----
df = merge(readsload, metadata, by="row.names", all = T)
## 将NA替换为0
df$readssum[is.na(df$readssum)] = 0
colnames(df)[1] = "Sample"

## 去除异常板的数据
# df = dplyr::filter(df, plate !=  opts$plate)

# 提取孔的行号，用于着色
split1 = as.data.frame(str_split_fixed(df$well, "",2))
colnames(split1)[1] = "row"
colnames(split1)[2] = "num"
df_split <- cbind(df, split1)

sample_A12 = dplyr::filter(df_split, well == opts$negative)
sample_B12 = dplyr::filter(df_split, well == opts$positive)


## 计算指定分位数的值
a = quantile(sample_A12$readssum,opts$threshold)

print(paste(opts$threshold*100,"%的阴性对照孔reads数低于：",a, sep=""))
print(paste("Reads of ", opts$threshold*100, "% negative control less than：",a, sep=""))

write.table(a, opts$output, sep = "\t", quote = F, col.names = F,row.names = F, na="")

## 绘制阴性对照箱线图#----
# p = ggplot(sample_A12, aes(well, readssum)) + 
#     geom_boxplot() + 
#     geom_hline(yintercept=a, linetype=2, color="blue") + 
#     annotate("text", x=0.6, y=a+0.05*max(sample_A12$readssum), label=a)
# ggsave(paste(opts$output, ".pdf", sep=""), p, width = 2, height = 3)
# ggsave(paste(opts$output, ".png", sep=""), p, width = 2, height = 3)

## 绘制阴、阳性对照箱线图#----
control = rbind(sample_A12, sample_B12)
control$well = factor(as.character(control$well), levels = unique(as.character(control$well)))
levels(control$well) = c("Negative", "Positive")
p = ggplot(control, aes(well, readssum, fill = well)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter() + 
  geom_hline(yintercept=a, linetype=2, color="blue") + 
  annotate("text", x=0.6, y=a+0.02*max(control$readssum), label=a, size = 2.5) + 
  main_theme + xlab("Controls")+ylab("Count of reads") + 
  theme(legend.position = "NA") # + coord_flip()
ggsave(paste(opts$output, ".control.pdf", sep=""), p, width = 89, height = 59, unit = "mm")
ggsave(paste(opts$output, ".control.png", sep=""), p, width = 89, height = 59, unit = "mm")

# 对数转换#----
# p = ggplot(control, aes(well, log2(readssum+1), fill = well)) + 
#   geom_boxplot(outlier.shape = NA) + geom_jitter() +
#   geom_hline(yintercept=log2(a+1), linetype=2, color="blue") + 
#   annotate("text", x=0.6, y=log2(a+1), label=a) 
# ggsave(paste(opts$output, ".control_log2.pdf", sep=""), p, width = 4, height = 3)
# ggsave(paste(opts$output, ".control_log2.png", sep=""), p, width = 4, height = 3)

print(paste("Result figure in ",opts$output, ".control.pdf/png", sep=""))

print("Analysis completed!!!")
