[TOC]

# 培养组学分析流程

作者：刘永鑫 (yxliu@genetics.ac.cn)等 (遗传发育所白洋组)

版本： v1.1   日期：2022-03-28

使用此流程，请引用:Jingying Zhang, Yong-Xin Liu, Xiaoxuan Guo, Yuan Qin, Ruben Garrido-Oter, Paul Schulze-Lefert & Yang Bai. 2021. High-throughput cultivation and identification of bacteria from the plant root microbiota. **Nature Protocols** 16: 988-1012. https://doi.org/10.1038/s41596-020-00444-7

使用方法：请先按照README.md中步骤安装culturome的conda环境，并下载Culturome项目。然后在RStudio中逐行运行本脚本，或复制缩进的代码在Terminal中运行。


## 分析步骤 

设置工作目录 ( work directory，`wd`)；修改 `db` 为 `Culturome`的安装目录
    
    # 启动culturome的conda环境
    conda activate culturome
    # 修改为工作目录，无d盘可改为c盘
    wd=/mnt/d/culture
    # 创建项目文件夹并进入，-p允许建立多级目录
    mkdir -p $wd && cd $wd
    # 设置Culturome位置
    db=~/Culturome
    # 建议项目所需文件夹
    mkdir -p seq temp result

(可选)下载测试数据并解压。或准备测序数据于seq目录中
    
    # 下载示例数据，-c支持断点续传。下载中断可再次运行此命令尝试
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_f1.fq.gz -O seq/L1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_r2.fq.gz -O seq/L1_2.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_f1.fq.gz -O seq/L2_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_r2.fq.gz -O seq/L2_2.fq.gz
    gunzip seq/*.gz

### 1. 编写样本信息(mapping files)

为每个测序文库编写mapping files。每个孔为双端标签的16S扩增子测序。

    # 生成文库列表，
    ls seq/*.fq | cut -f 2 -d '/' | cut -f 1 -d '_' | uniq > result/library.txt
        
    # 循环为每个文库生成48x96行的mapping file 
    # 方法1. 第个文库样本信息参数相同，使用for循环批量生成，不同需手动修改结果文件
    for l in `cat result/library.txt`; do
      write_mapping_file.pl -i ${db}/script/barcodeF96.txt -b ${db}/script/barcodeR48.txt \
        -F AACMGGATTAGATACCCKG -R ACGTCATCCCCACCTTCC \
        -L ${l} -p 48 -v Nippobare -c Root -m TSB -B 1 -s Rice -d WildType \
        -o seq/${l}.txt; done
    
    # 方法2. 库中参数不同时准备准备doc/library.txt列表填写各库中信息
    # awk -v nar="$db" 'BEGIN{OFS=FS="\t"}{system("write_mapping_file.pl -i "nar"/script/barcodeF96.txt -b "nar"/script/barcodeR48.txt -F AACMGGATTAGATACCCKG -R ACGTCATCCCCACCTTCC -L "$1" -p 48 -v "$5" -c "$6" -m "$4" -B "$7" -s "$3" -d "$9" -o seq/"$1".txt");}' <(tail -n+2 doc/library.txt)

    # 合并mapping file(s)为metadtata
    l=`tail -n+2 result/library.txt|cut -f1|head -n1`
    cat <(head -n1 seq/${l}.txt | sed 's/#//g') \
        <(cat seq/*.txt |grep -v '#'|grep -v -P '^SampleID\t') \
        > result/metadata.txt

### 2. 校验mapping file (可选)

手动编写或修改的文章必须校验，确保下游分析顺利完成

    for l in `cat result/library.txt`; do
      validate_mapping_file.py -m seq/${l}.txt -o temp/;done

### 3. 双端序列合并

    # time统计计算时间，rush指定并行任务数量，3s-3m15s
    time cat result/library.txt|rush -j 3 \
      'vsearch -fastq_mergepairs \
        seq/{1}_1.fq -reverse seq/{1}_2.fq \
	    -fastqout temp/{1}.fq'

### 4. 样本拆分

    # 提取barcodes, 1-15m, 注意根据实验设计修改barcode长度
    time cat result/library.txt|rush -j 3 \
      'extract_barcodes.py \
        -f temp/{1}.fq -m seq/{1}.txt \
    	  -c barcode_paired_stitched \
    	  --bc1_len 10 --bc2_len 6 \
    	  -a --rev_comp_bc2 \
    	  -o temp/{1}'

    # 样本拆分, 2m8s
    time cat result/library.txt|rush -j 3 \
      'split_libraries_fastq.py \
        -i temp/{1}/reads.fastq \
        -b temp/{1}/barcodes.fastq \
    	  -m seq/{1}.txt \
    	  -q 19 --max_barcode_errors 0 \
    	  --barcode_type 16 --phred_offset 33 \
    	  -o temp/{1}'

	# 格式化为vsearch输入文件
	rm -rf temp/qc.fa
  for l in `cat result/library.txt`; do
    cut -f 1 -d ' ' temp/${l}/seqs.fna | sed 's/_/./' \
	    >> temp/qc.fa; done

### 5. 统计每个板和孔中的测序读长

Visualize counts of samples in library. We using home-made script to visualize the reads distribution in bar plot.

    mkdir -p result/split
    for l in `cat result/library.txt`; do
    # 准备绘图文件
    tail -n+16 temp/${l}/split_library_log.txt| \
        head -n-4 > result/split/${l}.txt
    # 绘制每个文库的读长分布柱状图
    stat_split_bar.R \
        -i result/metadata.txt \
        -d result/split/${l}.txt \
        -o result/split/
    done

- "每个孔的测序数量 `result/split/L1.txt.well.pdf/png`"
- "每个孔的测序数量直方图 `result/split/L1.txt.histogram.pdf/png`"
- "每个板的测序数量柱状图 `result/split/L1.txt.plate.pdf/png`"

### 6. 切除引物

    # 1s-2m，根据实验设置引物长度
    time vsearch --fastx_filter temp/qc.fa \
      --fastq_stripleft 19 --fastq_stripright 18 \
      --fastaout temp/filtered.fa
      
### 7. 挑选代表序列

    # 计算每个读长出现的频率，去除小于10的读长， 4s-3m
    time vsearch --derep_fulllength temp/filtered.fa \
    	--relabel Uni --minuniquesize 10 --sizeout \
    	--output temp/uniques.fa
    # 查看非冗余序列结果大小
    ls -hs temp/uniques.fa

    # 采用unoise3算法去噪鉴定ASV
    vsearch --cluster_unoise temp/uniques.fa --minsize 10 --centroids temp/Zotus.fa
    # 从头去除嵌合体
    vsearch --uchime3_denovo temp/Zotus.fa --relabel ASV_ --nonchimeras result/ASV.fa
    

### 8. 构建ASV表

	# 比对读长到ASV进行定量，98.84% matched, 38s-30m
	time vsearch --usearch_global temp/filtered.fa \
	    --db result/ASV.fa \
        --otutabout temp/ASV_table.txt \
        --id 0.97

### 9. 假阳性率控制

扩增子测序容易受低丰度DNA污染，采用阴性对照控制假阳性结果。

    # 计算样阴、阳性对照的数据量，确定过滤假阳性的阈值
    negative_threshold.R \
        --input temp/ASV_table.txt --metadata result/metadata.txt \
        --threshold 1 \
        --negative A12 --positive B12 \
        --output result/fdr.txt
    cat result/fdr.txt

    # 过滤潜在假阳性孔
    # e.g. Deleted 705，Keep 4990
    otutab_trim.R \
        --input temp/ASV_table.txt \
        --min_sample_size `cat result/fdr.txt` \
        --output result/ASV_table.txt
    
    # 统计每个板中阳性孔数量
    head -n1 result/ASV_table.txt | cut -f2- | \
        sed 's/\t/\n/g' | cut -c1-5 | sort | \
        uniq -c | sort -k1,1n

### 10. 物种注释

    # 基于RDP训练集16，置信度阈值0.6.
    vsearch --sintax result/ASV.fa --db ${db}/db/rdp_16s_v16_sp.fa --tabbedout temp/ASV.fa.tax --sintax_cutoff 0.6
    # 制作2列，8列的物种注释表
    cut -f 1,4 temp/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy_2.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' |sed 's/#//g;s/ //g' > result/taxonomy_8.txt
    
    # ASV表按属合并
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$7} NR>FNR{print $0,a[$1]}' result/taxonomy_8.txt result/ASV_table.txt | sed '/\t$/d' | sed '1 s/Genus/KO/' > result/ASV_table7.txt
    mat_gene2ko.R -i result/ASV_table7.txt -o result/genus -n 1000000

### 11. 鉴定ASV非冗余的菌

    time identify_isolate.R \
        --input result/ASV_table.txt \
        --genus result/genus.count \
        --taxonomy result/taxonomy_8.txt \
        --output result/isolate

- "稀释曲线`result/isolate_rare_curve.pdf/png`"
- "孔列表`result/isolate_well.txt`"
- "ASV 列表`result/isolate_ASV.txt`"
- "ASV/G属纯度分布`result/isolate.Distribution.pdf/png`"

    # (可选) 每个板的纯度热图
    mkdir -p result/purity
    for l in `cat result/library.txt`; do
    # Prepare input file
    awk 'NR==FNR{a[$1]=$3}NR>FNR{print $1"\t"a[$1]}' result/isolate_well.txt seq/${l}.txt | \
        tail -n+2 | sed 's/\t$/\t0/' \
        > result/purity/${l}.txt
    # Format list into plate format
    format_list2plate.pl -i result/purity/${l}.txt \
        -o result/purity/${l}/
    # Batch plot each plate
    list=`ls result/purity/${l}/|cut -f 1 -d '.'|cut -f 2 -d 'P'|sort|uniq`
    for plate in $list;do 
        plot_pheatmap.sh -i result/purity/${l}/${l}P${plate}.plate \
            -o result/purity/${l}/${l}P${plate}.png
    done; done

- 每个板中每个孔纯度热图`result/purity/L*/*.png`

### 12. 培养菌进化分支树

    # 准备graphlan输入文件
    graphlan_prepare.R --input result/ASV_table.txt \
        --taxonomy result/taxonomy_8.txt \
        --abundance 0 --number 150 \
        --output result/graphlan/ 
    # 绘图
    graphlan_plot.sh -i $db -o result/graphlan

结果树图：result/graphlan/graphlan.pdf

### 13. (可选)项目清理

项目结束后，删除临时文件，节省空间

    # 压缩原始数据
    gzip seq/*.fq
    # 删除临时文件
    rm -r temp
    # 结果打包
    zip result.zip -r result/


