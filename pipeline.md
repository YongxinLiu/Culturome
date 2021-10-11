# Culturome pipeline v1.0

Authors: Yong-Xin Liu (yxliu@genetics.ac.cn), Yuan Qin (yqin@genetics.ac.cn)

Date: 2020-08-12

This is an pipeline for one sequencing library. For multiple sequencing libraries, please follow `pipeline_multiple.md`.

If use it, please cite: Jingying Zhang, [Yong-Xin Liu](http://bailab.genetics.ac.cn/YongxinLiuEn.html), Xiaoxuan Guo, Yuan Qin, Ruben Garrido-Oter, Paul Schulze-Lefert & Yang Bai. (2021). High-throughput cultivation and identification of bacteria from the plant root microbiota. ***Nature Protocols***, doi: https://doi.org/10.1038/s41596-020-00444-7

## Procedure

Modify `wd` to absolute directory of `Culturome`, then run the following script to initialize your environment. set your library ID in the variable `l`.

    # e.g. pipeline in my github directory
    wd=~/github/Culturome
    # e.g. pipeline in Virtualbox image
    wd=/home/qiime/Culturome-master
    PATH=${wd}/script:$PATH
    mkdir -p temp result
    l=L1

Download the demo data in the working directory. 

    # Download test data and database
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127979/CRR127979_f1.fq.gz -O L1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127979/CRR127979_r2.fq.gz -O L1_2.fq.gz
    wget -c http://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz
    gunzip *.gz

Bioinformatic analysis to identify cultivated bacteria (Fig. 4) ● Timing 1–3 days (~1 hour for 1 GB data). The exact time depends on the size of the dataset and computational experience. 

### 1. Prepare the mapping file to analyze the amplicon data.

The mapping file contains metadata such as the sample ID, forward and reverse barcodes, and forward and reverse primers. Although the mapping file can be manually generated according to the [QIIME mapping file format,](https://qiime.org/documentation/file_formats.html#metadata-mapping-files)it is tedious and may generate errors when determining the positions of wells and each plate using combined well and plate barcodes, given the large number of rows (4,320 for 45 cell culture plates). Here, we provide an in-house script (write_mapping_file.pl script) to generate a mapping file using input including the well and plate barcodes (-i, -b), forward and reverse primers (-F, -R), and other detailed information (-L, -p, -v, -c, -m, -s) (Supplementary Table 1 and 2). For detailed descriptions of the parameters, please type `write_mapping_file.pl -h`.

    write_mapping_file.pl \
        -i ${wd}/script/barcodeF96.txt \
        -b ${wd}/script/barcodeR48.txt \
        -F AACMGGATTAGATACCCKG \
        -R ACGTCATCCCCACCTTCC \
        -L ${l} -p 45 -v Nippobare \
        -c Root -m TSB -s Rice \
        -o ${l}.txt

![image](http://210.75.224.110/github/Culturome/script/fig/L1.jpg)

Table 1. Example of the mapping file. It must start with #SampleID. SampleID `L1P01A1` represent library 1 (L1), plate 1 (P1), and well A1. BarcodeSequence is forward plus reverse barcodes.

### 2. (Optional) Validate the mapping file 

There are many format requirements for mapping file, and errors will affect the following analysis. Using `validate_mapping_file.py` form QIIME to determine whether the mapping file is formatted properly.

    validate_mapping_file.py \
        -m ${l}.txt \
        -o temp/

 If the message "`No errors or warnings were found in mapping file.`" appear, your file is correct. Otherwise, revise your mapping file based on the output report (`*.html, *.log and *_corrected.txt`).
 
### 3. Merge paired-end reads into single-end reads

The 250-bp paired reads from HiSeq2500/NovaSeq6000 are merged according to the reverse complement of the end part of the sequences. 

    # vsearch merge pair-end reads, 5s, 98.9%
    time vsearch -fastq_mergepairs ${l}_1.fq \
        -reverse ${l}_2.fq \
	    -fastqout temp/${l}.fq 

### 4. Demultiplex barcoded sequences into samples corresponding to each well

According to the mapping file, remove plate and well barcodes (`extract_barcodes.py`), relabel sequences with the ID of plates and wells (`split_libraries_fastq.py`), and adjust the name of sequences to the format compatible with USEARCH. 

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

	# format to usearch compatible 
	cut -f 1 -d ' ' temp/${l}/seqs.fna | \
        sed 's/_/./' > temp/qc.fa

### 5. isualize sequence counts

Visualize sequence counts from each plate to check the success and evenness of amplicon sequencing. For more details, please type `stat_split_bar.R -h`.

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

**Fig. 1 | The bar plot showing read counts in each well of 96-well cell culture plates. n = 4320 (96  45).**


![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.histogram.png)

**Fig. 2 | The bar plot showing the distribution of read counts in each well.**

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.plate.png)

**Fig. 3 |  Bar plot showing the distribution of read counts in each plate generated by the two-sided barcode system**. Read counts of amplified 16S rRNA gene sequences in each plate of cultivated bacteria reveal the depth and evenness of Illumina sequencing. This figure includes sequencing data from 45 cell culture plates of cultivated bacteria derived from Nippornbare rice roots. The original data have been used in Nature biotechnology.

### 6. Remove forward and reverse primer sequences

    usearch -fastx_truncate temp/qc.fa \
        -stripleft 19 -stripright 18 \
        -fastaout temp/filtered.fa 

### 7. Identify amplicon sequence variants (ASVs)

Identify amplicon sequence variants (ASVs) by performing dereplication in VSEARCH and denoising in USEARCH. Final, rename Zotus to ASV.

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

Identify bacteria in each well of a 96-well cell culture plate by quantifying the ASV abundance by mapping all clean sequences to ASVs using VSEARCH.

	# 99% matched, 1~3m
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

Amplicon sequencing from each well is easily contaminated by low-abundant DNAs from reagents, the environment, or Illumina sequencing errors. Remove samples from wells containing fewer reads than the maximum number of reads in the negative controls. In each plate, A12 is the negative control with nuclease-free water and B12 is the positive control with E. coli DNA.

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

**Fig. 4 | Read counts for the negative (nuclease-free water as the PCR template) and positive controls (E. coli DNA as the PCR template).** The horizontal bars within boxes represent medians. The tops and bottoms of boxes represent the 75th and 25th percentiles, respectively. The upper and lower whiskers extend to data no more than 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. The dots represent the samples (*n* = 45).


### 10. Taxonomic classification

Classify the taxonomy of ASVs from cultivated bacteria based on the RDP train set 16 database and adjust to the tab-separated tabl

    # 30s, cutoff set 0.6
    usearch -sintax result/ASV.fa \
        -db rdp_16s_v16_sp.fa \
        -tabbedout temp/ASV.fa.tax \
        -sintax_cutoff 0.6 -strand both
    # summary phylum and genus, format to table
    tax_sum.sh -i temp/ASV.fa.tax \
        -d result/ASV_table.txt \
        -o result/

### 11. Identify non-redundancy isolates

Combining the ASV table and taxonomy, evaluate the saturation of bacterial ASV diversity based on the number of wells containing bacteria, obtain an overview of the distribution of wells containing different numbers of ASVs or genera, and examine the purity of the cultivated bacteria in each well. The outputs include two tables: an ASV list (`isolate_ASV.txt`, Supplementary Table 5) including five best wells containing bacteria with the corresponding 16S rRNA gene sequence and a well list (`isolate_well.txt`) including all detected ASV purity, read counts and taxonomy.


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

**Fig. 5 | Rarefaction curve of ASVs based on the number of wells containing bacteria.** The curve reaches the plateau stage, indicating that a sufficient number of plates were used for high-throughput bacterial isolation. The iteration number for each well number is 30. The horizontal bars within boxes represent medians. The tops and bottoms of boxes represent the 75th and 25th percentiles, respectively. The upper and lower whiskers extend to data no more than 1.5 × the interquartile range from the upper edge and lower edge of the box.

![image](http://210.75.224.110/github/Culturome/script/fig/isolate.DistributionPurityPer.png)

**Fig. 6 | Distribution of wells regarding the purity based on evaluation of ASVs of the bacterial 16S rRNA gene.** Consistent with the assumption that low abundance sequences may come from multiple polymorphic copies of the 16S rRNA gene and sequencing errors, we found that cultivated bacterial with a purity greater than 95% are likely to be pure when we performed streaking on agar plates. The ratio of wells with purify greater than 95% to the rest wells containing multiple bacteria is 3.6:1. 

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_ASV.jpg)

Table 2. Candidate wells and taxonomy of each ASV in ASV list (`isolate_ASV.txt`). Each ASV contained 5 top optimal candidates (Sorted by purity and abundance).

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_sample.jpg)

Table 3. Counts, purity and taxonomy in each well (`isolate_sample.txt`)

![image](http://210.75.224.110/github/Culturome/script/fig/purity.png)

**Figure 7. An example of the purity of wells containing cultivated bacteria in a 96-well cell culture plate.**


### 12. Summarize the taxonomy

Summarize the taxonomic distribution and occurrence frequency of cultivated bacteria using GraPhlAn. Default labeled in families. More detail type `graphlan_prepare.R -h`. If too much families leads to overlapping text, plese using `graphlan_prepare_order.R` instead of `graphlan_prepare.R` to label in order level.


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

**Fig 8. Cladogram showing the taxonomic distribution and occurrence frequency of cultivated bacteria.** The inner ring represents the dereplicated ASVs from cultivated root bacteria. The heat map in the outer ring represents the log2-transformed number of cultivated bacterial isolates belonging to the corresponding ASV.


### 13.	(Optional) Cross-reference the cultivated bacteria with 16S profile

Cross-reference the cultivated bacteria with the corresponding root microbiota profiling data using the similarity of V5–V7 regions in the 16S rRNA gene. 

    cross_reference.sh \
        -i script/profiling_ASV.fa \
        -d script/profiling_ASVtab.txt \
        -r result/ASV.fa \
        -o result/cross_reference.txt


### 14. Project cleanup

When the project is complete, you can compress input files and delete temporary files to save space. All the important results in `result` folder.

    # Compress raw data and database
    gzip *.fq *.fa
    # Delete temporary files
    rm -r temp
    # Packaging results
    zip script/result_single.zip -r result/

