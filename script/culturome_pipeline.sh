#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
# 程序目录
wd=/mnt/bai/yongxin/github/Culturome/
# 输入文件
input=/mnt/bai/yongxin/github/Culturome/db/demo_input.txt
# 输出目录绝对路径
output=/var/www/html/culturome/tmp/result/
# 输出目录名称
name=email+time

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    culturome_pipeline.sh
Version:     1.0
Date:        2020/4/12
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Description: Based on download link to analyze culturome data
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, 
Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, 
Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. 
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
Nature Biotechnology. 2019, 37: 676-684. 
doi:10.1038/s41587-019-0104-4
-------------------------------------------------------------------------------
Version 1.0 2020/4/12

OPTIONS:
	-i input download links, default ${input}
	-d directory of Culturome, default ${wd}
	-n name of temp project, default ${name}
	-o output directory, default ${output}
	-m metadatafile, default ${metadata}
	-? show help of script 显示帮助

Example:
culturome_pipeline.sh -i ${input} -d ${wd} -o ${output} -n ${name} -m ${metadata}

EOF
}


# 参数解析 Analysis parameter
while getopts "i:o:d:n:m:" OPTION
do
	case $OPTION in
		i)
			input=$OPTARG
			;;
		d)
			wd=$OPTARG
			;;
		n)
			name=$OPTARG
			;;
		m)
			metadata=$OPTARG
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
cfg=${wd}/cfg/
# 输入文件
seq=tmp/${name}/seq
# 临时文件
tmp=tmp/${name}/tmp

echo 'Please check the following parameters!'

echo ${metadata}
echo result/${name}


# 清理工作环境
rm -rf seq result
mkdir -p seq result

# 下载文件
# 基于列表批量下载并改名
awk 'BEGIN{OFS="\t";FS="[ \t]+"}{system("wget -c "$1" -O seq/"$2" ")}' ${input}
gunzip seq/*.gz


# 1. Generate library.txt from seq, fastq file in unzipped and extended in .fq
ls seq/*.fq|cut -f 2 -d '/'|cut -f 1 -d '_'|uniq \
    >result/library.txt
# Loop write mapping file. For different paramters of each library need write manual
for l in `cat result/library.txt`; do
write_mapping_file.pl \
    -i ${wd}/script/barcodeF96.txt \
    -b ${wd}/script/barcodeR48.txt \
    -F AACMGGATTAGATACCCKG \
    -R ACGTCATCCCCACCTTCC \
    -L ${l} -p 48 -v Nippobare \
    -c Root -m TSB -B 1 -s Rice -d WildType \
    -o seq/${l}.txt
done

# Merge mapping file(s) into one metadata
cat <(head -n1 seq/${l}.txt | sed 's/#//g') \
    <(cat seq/*.txt |grep -v '#'|grep -v -P '^SampleID\t') \
    > result/metadata.txt

# 外源mapping file，有则替换metadata.txt
if [ -f "${metadata}" ]; then
    sed -i 's/#//' ${metadata}
    mv ${metadata} result/metadata.txt
fi

# 2. (Optional)Validate mapping file


# 3. Merge pair-end reads, 5s
# mkdir -p tmp
for l in `cat result/library.txt`; do
vsearch -fastq_mergepairs seq/${l}_1.fq \
    -reverse seq/${l}_2.fq \
    -fastqout tmp/${l}.fq
done


# 4. Demultiplexing
# extract barcodes, 1m37s, note barcode legnth
rm -rf tmp/qc.fa
for l in `cat result/library.txt`; do
extract_barcodes.py \
    -f tmp/${l}.fq -m seq/${l}.txt \
	-c barcode_paired_stitched \
	--bc1_len 10 --bc2_len 6 \
	-a --rev_comp_bc2 \
	-o tmp/${l}
# split library, 2m8s
split_libraries_fastq.py \
    -i tmp/${l}/reads.fastq \
    -b tmp/${l}/barcodes.fastq \
    -m seq/${l}.txt \
    -q 19 --max_barcode_errors 0 \
    --barcode_type 16 --phred_offset 33 \
    -o tmp/${l}
# format to usearch
cut -f 1 -d ' ' tmp/${l}/seqs.fna \
    | sed 's/_/./' \
    >> tmp/qc.fa
done


# 5. Visualize counts of samples
# format split library
mkdir -p result/split
for l in `cat result/library.txt`; do
tail -n+16 tmp/${l}/split_library_log.txt| \
    head -n-4 > result/split/${l}.txt
# plot each library bar
stat_split_bar.R \
    -i result/metadata.txt \
    -d result/split/${l}.txt \
    -o result/split/
done


# 6. Cut primers
usearch -fastx_truncate tmp/qc.fa \
    -stripleft 19 -stripright 18 \
    -fastaout tmp/filtered.fa


# 7. Pick representitve sequences
# Calculate frequency of non-reduncancy reads, 4s
vsearch \
    --derep_fulllength tmp/filtered.fa \
    --relabel Uni --minuniquesize 8 --sizeout \
    --output tmp/uniques.fa 
# Denoise by unoise3, 2s
usearch -unoise3 tmp/uniques.fa \
    -zotus tmp/Zotus.fa
# Rename to ASV
awk 'BEGIN {n=1}; />/ {print ">ASV_" n; n++} !/>/ {print}' tmp/Zotus.fa \
    > result/ASV.fa


# 8. (Optional) Remove host and unspecific amplification


# 9. Construct ASV table
vsearch --usearch_global tmp/filtered.fa \
    --db result/ASV.fa \
    --otutabout tmp/ASV_table.txt \
    --id 0.97


# 10. Calculate 100%(1) false discovery reads, mean cut all negative control
negative_threshold.R \
    --input tmp/ASV_table.txt \
    --metadata result/metadata.txt \
    --threshold 1 \
    --negative A12 \
    --positive B12 \
    --output result/fdr.txt

# Filter flase discovery well in feature table
usearch -otutab_trim tmp/ASV_table.txt \
    -output result/ASV_table.txt \
    -min_sample_size `cat result/fdr.txt` 


# 11. Taxonomic classification
usearch -sintax result/ASV.fa \
    -db ${wd}/db/rdp_16s_v16_sp.fa \
	-tabbedout tmp/ASV.fa.tax \
	-sintax_cutoff 0.6 -strand both
tax_sum.sh -i tmp/ASV.fa.tax \
    -d result/ASV_table.txt \
    -o result/


## Plotting taxonomy(no response)
## Prepare graphlan files
#graphlan_prepare_order.R \
#    --input result/ASV_table.txt \
#    --output result/graphlan/ \
#    --taxonomy result/taxonomy_8.txt \
#    --abundance 0 \
#    --number 1000
## Plot graphlan
#conda deactivate
#graphlan_plot.sh -i ${wd} \
#    -o result/graphlan


# 12. Identify non-redundancy isolates
identify_isolate.R \
    --input result/ASV_table.txt \
    --taxonomy result/taxonomy_8.txt \
    --output result/isolate


# 13. Each well purity in plate
# Fetch purity in result table
mkdir -p result/purity
for l in `cat result/library.txt`; do
awk 'NR==FNR{a[$1]=$3}NR>FNR{print $1"\t"a[$1]}' result/isolate_well.txt seq/${l}.txt | tail -n+2 | sed 's/\t$/\t0/' \
    > result/purity/${l}.txt
# Format list into plate format
format_list2plate.pl -i result/purity/${l}.txt \
    -o result/purity/${l}/

# Batch plot each plate
list=`ls result/purity/${l}/|cut -f 1 -d '.'|cut -f 2 -d 'P'|sort|uniq|head -n3`
for plate in $list;do 
    plot_pheatmap.sh -i result/purity/${l}/${l}P${plate}.plate \
        -o result/purity/${l}/${l}P${plate}.png
done
done
rm -rf plot_pheatmap.r Rplots.pdf

# 输入文件
cat ${input} > result/input_file.txt
# 使用说明
cat /mnt/bai/yongxin/github/Culturome/README.md /mnt/bai/yongxin/github/Culturome/pipeline.md > result/Readme.md
# 打包结果
zip download/${name}.zip -r result/

# touch download/${name}.zip

# 清理项目中间文件
# rm -rf seq/ result/
