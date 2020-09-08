# Culturome pipeline v1.0 (multiple)

Authors: Yong-Xin Liu (yxliu@genetics.ac.cn), Yuan Qin (yqin@genetics.ac.cn)

Date: 2020-08-12

This is an example pipeline for multiple libraries.

If you use this pipeline, please cite:

Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. **Nature Biotechnology**. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4


## Procedure

Modify `wd` to absolute directory of `Culturome`, then run the following script to initial your environment.

    wd=/mnt/bai/yongxin/github/Culturome
    PATH=${wd}/script:$PATH
    rm -rf temp result
    mkdir -p seq temp result

Make sure you have downloaded `Culturome` directory, and set the right path. Input fastq (unzipped) files in `seq` folder, and database in current folder is ready. The following codes to download example data, and database.

    # Download two example libraries
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_f1.fq.gz -O seq/L1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_r2.fq.gz -O seq/L1_2.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_f1.fq.gz -O seq/L2_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_r2.fq.gz -O seq/L2_2.fq.gz
    gunzip seq/*.gz
    # Downlad database
    wget -c http://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz
    gunzip *.gz

### 1. Write the metadata mapping files

Write the metadata mapping files for each library. Each well is a standard 16S rDNA amplicon sequencing. Mapping file is the metadata of each well, including sample ID, forward and reverse barcodes, forward and reverse primers, species, date, location, and so on. You can manually write the mapping file according to [QIIME metadata mapping file format](https://qiime.org/documentation/file_formats.html#metadata-mapping-files). We recommend using `write_mapping_file.pl` script to generate mapping file, which depend on barcodes list files (`barcodeF96.txt` and `barcodeR48.txt`).

For multiple libraries, we need a library list file "result/library.txt". Then using the loop for each library

    # Generate library.txt from seq, fastq file in unzipped and extended in .fq
    ls seq/*.fq | cut -f 2 -d '/' | \
        cut -f 1 -d '_' | uniq \
        > result/library.txt
        
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
    
For detail description of the parameters, please type `write_mapping_file.pl -h`.

![image](http://210.75.224.110/github/Culturome/script/fig/L1.jpg)

Table 1. Example of the mapping file. Mapping file must start with #SampleID. SampleID `L1P01A1` represent library 1 (L1), plate 1 (P1), and well A1. BarcodeSequence is forward plus reverse barcodes.

### 2. (Optional) Validate the mapping file 

If you write the mapping files manually, validate the mapping file(s) format is requirement. There are many format requirements for mapping file, and errors will affect the following analysis. Using `validate_mapping_file.py` form QIIME to check whether the format of mapping file is OK. If show message "`No errors or warnings were found in mapping file.`" means your file is corrected. Otherwise to revise your mapping file according to output report (`*.html, *.log and *_corrected.txt`).

    for l in `cat result/library.txt`; do
    validate_mapping_file.py -m seq/${l}.txt \
        -o temp/
    done

### 3. Merge pair-end reads

Merge pair-end reads into single-end reads. Most amplicon sequencing using Illumina HiSeq2500/NovaSeq6000 platform on pair-end 250 bp mode. We first merged the pair-end into sing-end reads, according the complement of the reads end.

    for l in `cat result/library.txt`; do
    time vsearch -fastq_mergepairs seq/${l}_1.fq \
        -reverse seq/${l}_2.fq \
	    -fastqout temp/${l}.fq
	done

### 4. Demultiplexing

Demultiplexing means split library into samples. We sequenced 4608(48 plates * 96 wells) samples in one library, and pair-end barcodes to index each sample. We need to remove these barcodes and rename sequences according to mapping files. QIIME provides scripts to deal these problems. We can extract barcodes by `extract_barcodes.py`, and rename each sequences name according to Sample ID by `split_libraries_fastq.py`.

    for l in `cat result/library.txt`; do
    # extract barcodes, 1m37s, note barcode length
    extract_barcodes.py \
        -f temp/${l}.fq -m seq/${l}.txt \
    	-c barcode_paired_stitched \
    	--bc1_len 10 --bc2_len 6 \
    	-a --rev_comp_bc2 \
    	-o temp/${l}
    # split library, 2m8s
    split_libraries_fastq.py \
        -i temp/${l}/reads.fastq \
        -b temp/${l}/barcodes.fastq \
    	-m seq/${l}.txt \
    	-q 19 --max_barcode_errors 0 \
    	--barcode_type 16 --phred_offset 33 \
    	-o temp/${l}
	# format to usearch
	cut -f 1 -d ' ' temp/${l}/seqs.fna \
	    | sed 's/_/./' \
	    >> temp/qc.fa
	done

### 5. Summary counts of well and plate

Visualize counts of samples in library. We using home-made script to visualize the reads distribution in bar plot.

    # format split library
    mkdir -p result/split
    for l in `cat result/library.txt`; do
    tail -n+16 temp/${l}/split_library_log.txt| \
        head -n-4 > result/split/${l}.txt
    # plot each library bar
    stat_split_bar.R \
        -i result/metadata.txt \
        -d result/split/${l}.txt \
        -o result/split/
    done

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

Remove the forward and reverse primers by length.

    # length of primers used, 10s
    time usearch -fastx_truncate temp/qc.fa \
    	-stripleft 19 -stripright 18 \
    	-fastaout temp/filtered.fa 

### 7. Pick representative sequences

Pick representative sequences. We first remove the redundancy of all reads, and calculate the frequency of reads. Then using unoise3 algorithm to denoise into amplicon sequence variants (ASV).

    # Calculate frequency of non-redundancy reads, 4s
    vsearch \
        --derep_fulllength temp/filtered.fa \
    	--relabel Uni --minuniquesize 8 --sizeout \
    	--output temp/uniques.fa 

    # Denoise by unoise3, 2s
    usearch -unoise3 temp/uniques.fa \
        -zotus temp/Zotus.fa
    
    # Rename to ASV
    awk 'BEGIN {n=1}; />/ {print ">ASV_" n; n++} !/>/ {print}' temp/Zotus.fa \
        > result/ASV.fa

### 8. Construct ASV table

Construct ASV table. Finally, we map all clean amplicon against the ASV to quantify the frequency in each sample.
 
	# 99.07% matched, 38s
	vsearch --usearch_global temp/filtered.fa \
	    --db result/ASV.fa \
        --otutabout temp/ASV_table.txt \
        --id 0.97
    
    # Stat, e.g. 2678 samples, 378 OTUs
    usearch -otutab_stats temp/ASV_table.txt \
        -output temp/ASV_table.stat
    cat temp/ASV_table.stat

### 9. False discovery rate control

False discovery rate control. Amplicon sequencing of each well is easily contaminated by low abundant DNAs from reagents, air or Illumina sequencing errors. Remove samples from wells containing reads lower than the maximal reads in negative controls (Fig. 5b). In each plate, A12 is a negative control with sterile water, and B12 is a positive control with E. coli DNA.

    # Calculate 100%(1) false discovery reads, mean cut all negative control
    negative_threshold.R \
        --input temp/ASV_table.txt \
        --metadata result/metadata.txt \
        --threshold 1 \
        --negative A12 \
        --positive B12 \
        --output result/fdr.txt

    # Filter flase discovery well in feature table
    # e.g. Deleted 704 / 5696 samples
    usearch -otutab_trim temp/ASV_table.txt \
        -output result/ASV_table.txt \
        -min_sample_size `cat result/fdr.txt` 
    
    # Well number in each plate
    head -n1 result/ASV_table.txt | cut -f2- | \
        sed 's/\t/\n/g' | cut -c1-5 | sort | \
        uniq -c | sort -k1,1n

![image](http://210.75.224.110/github/Culturome/result/fdr.txt.control.png)

![image](http://210.75.224.110/github/Culturome/script/fig/fdr.txt.control.png)

**Fig. 4 | Boxplot showing the read counts in negative (sterile water as PCR template) and positive controls (E. coli DNA as PCR template).** The horizontal bars within boxes represent medians. The tops and bottoms of boxes represent the 75th and 25th percentiles, respectively. The upper and lower whiskers extend to data no more than 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. The dots represent the samples (*n* = 96).


### 10. Taxonomic classification

Taxonomic classification. Based on RDP train set 16 databases, we use sintax to classify taxonomy of ASV. The confidence cutoff set ot 0.6.

    usearch -sintax result/ASV.fa \
	    -db ${wd}/rdp_16s_v16_sp.fa \
    	-tabbedout temp/ASV.fa.tax \
    	-sintax_cutoff 0.6 -strand both
    # summary phylum and genus, format to table, 3m
    tax_sum.sh -i temp/ASV.fa.tax \
        -d result/ASV_table.txt \
        -o result/

### 11. Identify non-redundancy isolates.

Combining ASV table and taxonomy, evaluate the saturation of bacterial ASV diversity according to the number of wells containing bacteria (Fig. 5), overview the distribution of wells containing different numbers of ASVs or genera (Fig. 6), and examine the purity of cultivated bacteria in each well (Fig. 7). The outputs include two tables: an ASV list (isolate_ASV.txt) including five wells containing bacteria having corresponding 16S rRNA sequence; a well list (isolate_well.txt) including all detected ASV sequences and taxonomy. 

    # Claculate tables and figures, 2m
    identify_isolate.R \
        --input result/ASV_table.txt \
        --taxonomy result/taxonomy_8.txt \
        --output result/isolate

- "Rarefaction curve in `result/isolate_rare_curve.pdf/png`"
- "Well list in `result/isolate_well.txt`"
- "Top1 ASV number:" 316
- "ASV list in `result/isolate_ASV.txt`"
- "Distribution of ASV/Genus number of all wells in `result/isolate.Distribution.pdf/png`"

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_rare_curve.png)

**Fig. 5 | A rarefaction curve of ASVs according to the number of wells containing bacteria.** The curve reach plateau stage indicating that the number of plates in high-throughput bacterial isolation is sufficient. 

![image](http://210.75.224.110/github/Culturome/script/fig/isolate.Distribution.png)

**Fig. 6 | Distribution of ASV/Genus number of all wells.** Because pure bacteria may contain 16S rRNA gene variations resulting in different ASVs, wells containing more one ASVs belonging to single genus may be pure. 

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_ASV.jpg)

Table 2. Candidate wells and taxonomy of each ASV in ASV list (`isolate_ASV.txt`). Each ASV contained 5 top optimal candidates (Sorted by purity and abundance).

![image](http://210.75.224.110/github/Culturome/script/fig/isolate_sample.jpg)

Table 3. Counts, purity and taxonomy in each well (`isolate_sample.txt`)

    # Each well purity in plate
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

- Purify heatmap in `result/purity/L*/*.png`

![image](http://210.75.224.110/github/Culturome/script/fig/purity.png)

**Figure 7. An example of purity of wells containing cultivated bacteria in 96-well microtiter plate.**

### 12.Summarize the taxonomy

Summarize the taxonomic distribution and occurrence frequency of cultivated bacteria using GraphLan. Default labeled in families. More detail type `graphlan_prepare.R -h`. If too much families leads to overlapping text, plese using `graphlan_prepare_order.R` instead of `graphlan_prepare.R` to label in order level. `--number 150` threshold top 150 most abundance ASVs.

    # Prepare graphlan files
    graphlan_prepare.R \
        --input result/ASV_table.txt \
        --taxonomy result/taxonomy_8.txt \
        --output result/graphlan/ \
        --abundance 0 \
        --number 150
    # Plot graphlan
    graphlan_plot.sh -i ${wd} \
        -o result/graphlan

![image](http://210.75.224.110/github/Culturome/script/fig/graphlan.png)

**Fig 8. Cladogram showing the taxonomic distribution and occurrence frequency of cultivated bacteria.** The inner ring represents the dereplicated ASVs from cultivated root bacteria. Heat map in the outer ring represent the log2 transformed number of cultivated bacterial isolates belong to the corresponding ASV.

### 13. Project cleanup

When the project is complete, you can compress input files and delete temporary files to save space.

    # Compress raw data and database
    gzip seq/*.fq *.fa
    # Delete temporary files
    rm -r temp
    # Packaging results
    zip script/result_multiple.zip -r result/

All the important results in `result` folder.
