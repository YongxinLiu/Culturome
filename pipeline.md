[TOC]

# Culturome pipeline v1.0

Authors: Yong-Xin Liu (yxliu@genetics.ac.cn), Yuan Qin (yqin@genetics.ac.cn)

Date: 2020-04-21

This is an example pipeline for one library. As for multiple libraries, please followed `pipeline_multiple.md`.

If you use this pipeline, please cite:
Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

## Procedure

Modify `wd` to absolute directory of `Culturome`, then run the following script to initial your environment. `l` set your library ID.

    wd=/mnt/bai/yongxin/github/Culturome
    PATH=${wd}/script:$PATH
    mkdir -p temp result
    l=L1

Make sure you have downloaded `Culturome` directory, and set the right path. Input fastq (unzipped) file and database current folder is ready.

    # Download test data and database
    #wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127979/CRR127979_f1.fq.gz -O L1_1.fq.gz
    #wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127979/CRR127979_r2.fq.gz -O L1_2.fq.gz
    #wget -c http://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz
    #gunzip *.gz
    
### 1. Prepare the mapping file to analysis the amplicon data

Mapping file is the metadata containing sample ID, forward and reverse barcodes, forward and reverse primers, and so on. You can manually write the mapping file according to [QIIME metadata mapping file format](https://qiime.org/documentation/file_formats.html#metadata-mapping-files). Here, we provide a self-written script (write_mapping_file.pl script) to generate mapping file, using the input of the well and plate barcodes (-i, -b), forward and reverse primers (-F, -R), and other detailed information (-L, -p, -v, -c, -m, -s). For detailed descriptions of parameters, please type `write_mapping_file.pl -h`.

    write_mapping_file.pl \
        -i ${wd}/script/barcodeF96.txt \
        -b ${wd}/script/barcodeR48.txt \
        -F AACMGGATTAGATACCCKG \
        -R ACGTCATCCCCACCTTCC \
        -L ${l} -p 45 -v Nippobare \
        -c Root -m TSB -s Rice \
        -o ${l}.txt

![image](http://210.75.224.110/github/Culturome/script/fig/L1.jpg)

Table 1. Example of the mapping file. Mapping file must start with #SampleID. SampleID `L1P01A1` represent library 1 (L1), plate 1 (P1), and well A1. BarcodeSequence is forward plus reverse barcodes.

### 2. (Optional) Validate the mapping file 

There are many format requirements for mapping file, and errors will affect the following analysis. Using `validate_mapping_file.py` form QIIME to check whether the format of mapping file is OK. If show message "`No errors or warnings were found in mapping file.`" means your file is corrected. Otherwise to revise your mapping file according to output report (`*.html, *.log and *_corrected.txt`).

    validate_mapping_file.py \
        -m ${l}.txt \
        -o temp/

### 3. Merge paired-end reads

Merge paired-end reads into single-end reads. Most amplicon sequencing using Illumina HiSeq2500/NovaSeq6000 platform in paired-end 250 bp mode. We first merged the paired-end into single-end reads, according the reverse complement of the end part of reads.

    # vsearch merge pair-end reads, 5s, 98.9%
    time vsearch -fastq_mergepairs ${l}_1.fq \
        -reverse ${l}_2.fq \
	    -fastqout temp/${l}.fq 

### 4. Demultiplexing

Demultiplex barcoded sequences into sequences of each well. According to the mapping file, remove plate and well barcodes (`extract_barcodes.py`), relabel sequences with the ID of plates and wells (`split_libraries_fastq.py`), and adjust the name of sequences to the format compatible with USEARCH. 

    # 1m37s, barcode legnth in mapping file
    extract_barcodes.py \
        -f temp/${l}.fq -m ${l}.txt \
    	-c barcode_paired_stitched \
    	--bc1_len 10 --bc2_len 6 \
    	-a --rev_comp_bc2 \
    	-o temp/${l}

    # 2m8s, relabel reads
    split_libraries_fastq.py \
        -i temp/${l}/reads.fastq \
        -b temp/${l}/barcodes.fastq \
    	-m ${l}.txt \
    	-q 19 --max_barcode_errors 0 \
    	--barcode_type 16 --phred_offset 33 \
    	-o temp/${l}

	# format to usearch
	cut -f 1 -d ' ' temp/${l}/seqs.fna | \
	    sed 's/_/./' > temp/qc.fa

### 5. Visualize counts

Visualize sequence counts from each plate to check success and evenness of amplicon sequencing (Fig. 1). For more detail, please type `stat_split_bar.R -h`.

    # When outputting multiple files, create subdirectories
    mkdir -p result/split
    # Adjust input file format
    tail -n+16 temp/${l}/split_library_log.txt| \
        head -n-4 > result/split/${l}.txt
    # Visualize counts of well and plate
    stat_split_bar.R \
        --input L1.txt \
        --database result/split/${l}.txt \
        --output result/split/
        
- "Well counts in `result/split/L1.txt.well.pdf/png`"
- "Histogram of well counts in `result/split/L1.txt.histogram.pdf/png`"
- "Plate counts in `result/split/L1.txt.plate.pdf/png`"

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.well.png)

**Fig. 1 | The bar plot showing the distribution of read counts in each well. n = 4609 (96 X 48).**


![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.histogram.png)

**Fig. 2 | The bar plot showing the distribution of read counts in each wells.**

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.plate.png)

**Fig. 3 | The bar plot showing the distribution of read counts generated by two-side barcode system in each plate.** Reads count of amplified 16S rRNA gene sequences in each plate of cultivated bacteria reveal the depth and evenness of Illumina sequencing. n = 45.

### 6. Remove primers

Remove forward and reverse primers of sequences, length according to primers.

    # 5-30s
    usearch -fastx_truncate temp/qc.fa \
    	-stripleft 19 -stripright 18 \
    	-fastaout temp/filtered.fa 

### 7. Identify amplicon sequence variants (ASVs)

Remove the redundancy of all reads, and calculate the frequency of reads. Then using unoise3 algorithm to denoise into amplicon sequence variants (ASVs). Final, rename ASV.

    # Calculate frequency of reads, 4s
    vsearch \
        --derep_fulllength temp/filtered.fa \
    	--relabel Uni --minuniquesize 8 --sizeout \
    	--output temp/uniques.fa 

    # Denoise by unoise3, 2s
    usearch -unoise3 temp/uniques.fa \
        -zotus temp/Zotus.fa
    
    # Rename to ASV
    format_ASVID.sh -i temp/Zotus.fa \
        -o result/ASV.fa

### 8. Construct ASV table

Identify bacteria in each well of 96-well microtiter plates, by quantifying the ASV abundance through mapping all clean sequences to ASVs using VSEARCH. 

	# 99% matched, 1m
	vsearch \
	    --usearch_global temp/filtered.fa \
	    --db result/ASV.fa \
        --otutabout temp/ASV_table.txt \
        --id 0.97
    
    # Stat, e.g. 2694 Samples, 378 ASVs
    usearch -otutab_stats temp/ASV_table.txt \
        -output temp/ASV_table.stat
    cat temp/ASV_table.stat

### 9. Remove wells only containing background noise

False discovery rate control. Amplicon sequencing of each well is easily contaminated by low abundant DNAs from reagents, air or Illumina sequencing errors. Remove samples from wells containing reads lower than the maximal reads in negative controls (Fig. 5b). In each plate, A12 is a negative control with sterile water, and B12 is a positive control with E. coli DNA.

    # Remove 100%(1) false discovery wells
    negative_threshold.R \
        --input temp/ASV_table.txt \
        --metadata L1.txt \
        --threshold 1 \
        --negative A12 \
        --positive B12 \
        --output result/fdr.txt

    # Filter flase discovery well in feature table
    # Deleted 621 / 2694 samples
    usearch -otutab_trim temp/ASV_table.txt \
        -output result/ASV_table.txt \
        -min_sample_size `cat result/fdr.txt` 

![image](http://210.75.224.110/github/Culturome/script/fig/fdr.txt.control.png)

**Fig. 4 | Boxplot showing the read counts in negative (sterile water as PCR template) and positive controls (E. coli DNA as PCR template).** The horizontal bars within boxes represent medians. The tops and bottoms of boxes represent the 75th and 25th percentiles, respectively. The upper and lower whiskers extend to data no more than 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. The dots represent the samples (*n* = 45).


### 10. Taxonomic classification

Based on RDP train set 16 databases, use sintax to classify taxonomy of ASV.

    # 30s, cutoff set 0.6
    usearch -sintax result/ASV.fa \
	    -db ${wd}/db/rdp_16s_v16_sp.fa \
    	-tabbedout temp/ASV.fa.tax \
    	-sintax_cutoff 0.6 -strand both
    # summary phylum and genus, format to table
    tax_sum.sh -i temp/ASV.fa.tax \
        -d result/ASV_table.txt \
        -o result/

### 11. Identify non-redundancy isolates

Combining ASV table and taxonomy, evaluate the saturation of bacterial ASV diversity according to the number of wells containing bacteria (Fig. 5), overview the distribution of wells containing different numbers of ASVs or genera (Fig. 6), and examine the purity of cultivated bacteria in each well (Fig. 7). The outputs include two tables: an ASV list (isolate_ASV.txt) including five wells containing bacteria having corresponding 16S rRNA sequence; a well list (isolate_well.txt) including all detected ASV sequences and taxonomy. 

    # Claculate tables and figures, 2m
    identify_isolate.R \
        --input result/ASV_table.txt \
        --taxonomy result/taxonomy_8.txt \
        --output result/isolate
    # plot purity heatmap
    plot_purity.sh \
        -i result/isolate_well.txt \
        -l L1 \
        -o result/purity

- "Rarefaction curve in `result/isolate_rare_curve.pdf/png`"
- "Well list in `result/isolate_well.txt`"
- "Top1 ASV number:" 316
- "ASV list in `result/isolate_ASV.txt`"
- "Distribution of ASV/Genus number of all wells in `result/isolate.Distribution.pdf/png`"
- Purify heatmap in `result/purity/L1/*.png`

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_rare_curve.png)

**Fig. 5 | A rarefaction curve of ASVs according to the number of wells containing bacteria.** The curve reach plateau stage indicating that the number of plates in high-throughput bacterial isolation is sufficient. 

![image](http://210.75.224.110/github/Culturome/script/fig/isolate.Distribution.png)

**Fig. 6 | Distribution of ASV/Genus number of all wells.** Because pure bacteria may contain 16S rRNA gene variations resulting in different ASVs, wells containing more one ASVs belonging to single genus may be pure. 

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_ASV.jpg)

Table 2. Candidate wells and taxonomy of each ASV in ASV list (`isolate_ASV.txt`). Each ASV contained 5 top optimal candidates (Sorted by purity and abundance).

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_sample.jpg)

Table 3. Counts, purity and taxonomy in each well (`isolate_sample.txt`)

![image](http://210.75.224.110/github/Culturome/script/fig/purity.png)

**Figure 7. An example of purity of wells containing cultivated bacteria in 96-well microtiter plate.**


### 12.Summarize the taxonomy

Summarize the taxonomic distribution and occurrence frequency of cultivated bacteria using GraphLan. Default labeled in families. More detail type `graphlan_prepare.R -h`. If too much families leads to overlapping text, plese using `graphlan_prepare_order.R` instead of `graphlan_prepare.R` to label in order level.

    # Prepare graphlan files
    graphlan_prepare_order.R \
        --input result/ASV_table.txt \
        --taxonomy result/taxonomy_8.txt \
        --output result/graphlan/ \
        --abundance 0 \
        --number 1000
    # Plot graphlan
    graphlan_plot.sh -i ${wd} \
        -o result/graphlan

![image](http://210.75.224.110/github/Culturome/script/fig/graphlan.png)

**Fig 8. Cladogram showing the taxonomic distribution and occurrence frequency of cultivated bacteria.** The inner ring represents the dereplicated ASVs from cultivated root bacteria. Heat map in the outer ring represent the log2 transformed number of cultivated bacterial isolates belong to the corresponding ASV.

### 13. Project cleanup

When the project is complete, you can compress input files and delete temporary files to save space.

    # Compress raw data and database
    gzip *.fq *.fa
    # Delete temporary files
    rm -r temp
    # Packaging results
    zip script/result_single.zip -r result/

All the important results in `result` folder.
