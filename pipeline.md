[TOC]

# Culturome pipeline

Authors: Yong-Xin Liu (yxliu@genetics.ac.cn) et al. from Bailab IGDB, CAS

Version: 2.0    Date: 2022-03-28

If you use this pipeline, please cite: Jingying Zhang, Yong-Xin Liu, Xiaoxuan Guo, Yuan Qin, Ruben Garrido-Oter, Paul Schulze-Lefert, et al. 2021. High-throughput cultivation and identification of bacteria from the plant root microbiota. **Nature Protocols** 16: 988-1012. https://doi.org/10.1038/s41596-020-00444-7

For install Culturome, please see th appendix 1 for Windows/Linux, and appendix 2 for Mac. The default directory for Windows Subsystem for Linux. For Linux or Mac, change to directory as your environment.

## Procedure 

Setup work directory (`wd`) to your project.

Modify `db` to absolute directory of `Culturome`, then run the following script to initial your environment.
    
    # Open culturome conda 
    conda activate culturome
    # set work directory, default for Windows, Mac change to ~/Documents/microbiome/culture
    wd=/mnt/c/culture
    mkdir -p $wd && cd $wd
    # set Culturome path, Mac change to ~/Documents/microbiome/Culturome
    db=/mnt/c/microbiome/Culturome
    mkdir -p seq temp result

(Option: test data)Make sure you have downloaded `Culturome` directory, and set the right path. Input fastq (unzipped) files in `seq` folder, and database in current folder is ready. The following codes to download example data, and database.

    # Download two example libraries
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_f1.fq.gz -O seq/L1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_r2.fq.gz -O seq/L1_2.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_f1.fq.gz -O seq/L2_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127982/CRR127982_r2.fq.gz -O seq/L2_2.fq.gz
    gunzip seq/*.gz

### 1. Write the metadata mapping files

Write the metadata mapping files for each library. Each well is a standard 16S rDNA amplicon sequencing. Mapping file is the metadata of each well, including sample ID, forward and reverse barcodes, forward and reverse primers, species, date, location, and so on. You can manually write the mapping file according to [QIIME metadata mapping file format](https://qiime.org/documentation/file_formats.html#metadata-mapping-files). We recommend using `write_mapping_file.pl` script to generate mapping file, which depend on barcodes list files (`barcodeF96.txt` and `barcodeR48.txt`).

For multiple libraries, we need a library list file "result/library.txt". Then using the loop for each library

    # Generate library.txt from seq, fastq file in unzipped and extended in .fq
    ls seq/*.fq | cut -f 2 -d '/' | \
        cut -f 1 -d '_' | uniq \
        > result/library.txt
        
    # Loop write mapping file. 
    # Method 1. For same paramters of each library need write manual
    for l in `cat result/library.txt`; do
    write_mapping_file.pl \
        -i ${db}/script/barcodeF96.txt \
        -b ${db}/script/barcodeR48.txt \
        -F AACMGGATTAGATACCCKG \
        -R ACGTCATCCCCACCTTCC \
        -L ${l} -p 48 -v Nippobare \
        -c Root -m TSB -B 1 -s Rice -d WildType \
        -o seq/${l}.txt
    done
    
    # Method 2. For different paramters of each library
    # awk -v nar="$db" 'BEGIN{OFS=FS="\t"}{system("write_mapping_file.pl -i "nar"/script/barcodeF96.txt -b "nar"/script/barcodeR48.txt -F AACMGGATTAGATACCCKG -R ACGTCATCCCCACCTTCC -L "$1" -p 48 -v "$5" -c "$6" -m "$4" -B "$7" -s "$3" -d "$9" -o seq/"$1".txt");}' <(tail -n+2 doc/library.txt)

    # Merge mapping file(s) into one metadata
    l=`tail -n+2 result/library.txt|cut -f1|head -n1`
    cat <(head -n1 seq/${l}.txt | sed 's/#//g') \
        <(cat seq/*.txt |grep -v '#'|grep -v -P '^SampleID\t') \
        > result/metadata.txt
    
For detail description of the parameters, please type `write_mapping_file.pl -h`.

![image](http://210.75.224.110/github/Culturome/script/fig/L1.jpg)

Table 1. Example of the mapping file. Mapping file must start with #SampleID. SampleID `L1P01A1` represent library 1 (L1), plate 1 (P1), and well A1. BarcodeSequence is forward plus reverse barcodes.

### 2. (Optional) Validate the mapping file 

If you write the mapping files manually, validate the mapping file(s) format is requirement. There are many format requirements for mapping file, and errors will affect the following analysis. Using `validate_mapping_file.py` form QIIME to check whether the format of mapping file is OK. If show message "`No errors or warnings were found in mapping file.`" means your file is corrected. Otherwise to revise your mapping file according to output report (`*.html, *.log and *_corrected.txt`).

    for l in `cat result/library.txt`; do
      validate_mapping_file.py -m seq/${l}.txt -o temp/;done

### 3. Merge pair-end reads

Merge pair-end reads into single-end reads. Most amplicon sequencing using Illumina HiSeq2500/NovaSeq6000 platform on pair-end 250 bp mode. We first merged the pair-end into sing-end reads, according the complement of the reads end. If error cmd, change '-j 3' to '-j 1'

    # Parallel, 3s, 3m15s
    time cat result/library.txt|rush -j 3 \
      'vsearch -fastq_mergepairs \
        seq/{1}_1.fq -reverse seq/{1}_2.fq \
	    -fastqout temp/{1}.fq'

### 4. Demultiplexing

Demultiplexing means split library into samples. We sequenced 4608(48 plates * 96 wells) samples in one library, and pair-end barcodes to index each sample. We need to remove these barcodes and rename sequences according to mapping files. QIIME provides scripts to deal these problems. We can extract barcodes by `extract_barcodes.py`, and rename each sequences name according to Sample ID by `split_libraries_fastq.py`.

    # extract barcodes, 1-15m, note barcode length
    time cat result/library.txt|rush -j 3 \
      'extract_barcodes.py \
        -f temp/{1}.fq -m seq/{1}.txt \
    	-c barcode_paired_stitched \
    	--bc1_len 10 --bc2_len 6 \
    	-a --rev_comp_bc2 \
    	-o temp/{1}'

    # split library, 2m8s
    time cat result/library.txt|rush -j 3 \
      'split_libraries_fastq.py \
        -i temp/{1}/reads.fastq \
        -b temp/{1}/barcodes.fastq \
    	-m seq/{1}.txt \
    	-q 19 --max_barcode_errors 0 \
    	--barcode_type 16 --phred_offset 33 \
    	-o temp/{1}'

	# format to usearch
	rm -rf temp/qc.fa
    for l in `cat result/library.txt`; do
	cut -f 1 -d ' ' temp/${l}/seqs.fna \
	    | sed 's/_/./' \
	    >> temp/qc.fa; done

### 5. Summary counts of well and plate

Visualize counts of samples in library. We using home-made script to visualize the reads distribution in bar plot.

    # plot each library bar
    mkdir -p result/split
    for l in `cat result/library.txt`; do
      tail -n+16 temp/${l}/split_library_log.txt| \
        head -n-4 > result/split/${l}.txt
      stat_split_bar.R \
        -i result/metadata.txt \
        -d result/split/${l}.txt \
        -o result/split/
    done

- Well counts in `result/split/L1.txt.well.pdf/png`
- Histogram of well counts in `result/split/L1.txt.histogram.pdf/png`"
- "Plate counts in `result/split/L1.txt.plate.pdf/png`"
- Each figure include PDF and PNG format, L1.txt and L1.txt.plate.txt are raw data for figures

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.histogram.png)

**Fig. 1 | The bar plot showing the distribution of read counts in each wells.**

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.plate.png)

**Fig. 2 | The bar plot showing the distribution of read counts generated by two-side barcode system in each plate.** Reads count of amplified 16S rRNA gene sequences in each plate of cultivated bacteria reveal the depth and evenness of Illumina sequencing. n = 45.

![image](http://210.75.224.110/github/Culturome/script/fig/L1.txt.well.png)

**Fig. 3 | The bar plot showing the distribution of read counts in each well. n = 4608 (96 X 48).**

### 6. Remove primers

Remove the forward and reverse primers by length.

    # 1s-2m
    time vsearch --fastx_filter temp/qc.fa \
      --fastq_stripleft 19 --fastq_stripright 18 \
      --fastaout temp/filtered.fa
      
### 7. Pick representative sequences

Pick representative sequences. We first remove the redundancy of all reads, and calculate the frequency of reads. Then using unoise3 algorithm to denoise into amplicon sequence variants (ASV).

    # 7.1 Calculate frequency of reads, 4s-3m
    time vsearch --derep_fulllength temp/filtered.fa \
    	--relabel Uni --minuniquesize 20 --sizeout \
    	--output temp/uniques.fa
    # check size, MB is OK, GB maybe wrong
    ls -hs temp/uniques.fa

    # 7.2 Denoise by unoise3 of vsearch, 2s
    vsearch --cluster_unoise temp/uniques.fa --minsize 20 --centroids temp/Zotus.fa
    
    # 7.3 De novo remove chimera 
    vsearch --uchime3_denovo temp/Zotus.fa --relabel ASV_ --nonchimeras result/ASV.fa
    

### 8. Construct ASV table

Construct ASV table. Finally, we map all clean amplicon against the ASV to quantify the frequency in each sample.
 
	# create feature table, 38s-30m
	time vsearch --usearch_global temp/filtered.fa \
	    --db result/ASV.fa \
        --otutabout temp/ASV_table.txt \
        --id 0.97

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
    cat result/fdr.txt

    # Filter flase discovery well in feature table
    # e.g. Deleted 691 / 5676; 8313/30601 samples
    otutab_trim.R \
        --input temp/ASV_table.txt \
        --min_sample_size `cat result/fdr.txt` \
        --output result/ASV_table.txt
    
    # Well number in each plate
    head -n1 result/ASV_table.txt | cut -f2- | \
        sed 's/\t/\n/g' | cut -c1-5 | sort | \
        uniq -c | sort -k1,1n > result/plate_positive.count
    cat result/plate_positive.count

![image](http://210.75.224.110/github/Culturome/script/fig/fdr.txt.control.png)

**Fig. 4 | Boxplot showing the read counts in negative (sterile water as PCR template) and positive controls (E. coli DNA as PCR template).** 


### 10. Taxonomic classification

Taxonomic classification. Based on RDP train set 16 databases, we use sintax to classify taxonomy of ASV. The confidence cutoff set ot 0.6.

    vsearch --sintax result/ASV.fa --db ${db}/db/rdp_16s_v16_sp.fa --tabbedout temp/ASV.fa.tax --sintax_cutoff 0.6
    
    # summary phylum and genus, format to table, 3-30m
    cut -f 1,4 temp/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy_2.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' |sed 's/#//g;s/ //g' > result/taxonomy_8.txt
    # ASV table to genus table
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$7} NR>FNR{print $0,a[$1]}' \
      result/taxonomy_8.txt result/ASV_table.txt | \
    	sed '/\t$/d' | sed '1 s/Genus/KO/' > result/ASV_table7.txt
    # 转换为基因表
    mat_gene2ko.R -i result/ASV_table7.txt -n 100 \
      -o result/genus

### 11. Identify non-redundancy isolates

Combining ASV table and taxonomy, evaluate the saturation of bacterial ASV diversity according to the number of wells containing bacteria (Fig. 5), overview the distribution of wells containing different numbers of ASVs or genera (Fig. 6), and examine the purity of cultivated bacteria in each well (Fig. 7). The outputs include two tables: an ASV list (isolate_ASV.txt) including five wells containing bacteria having corresponding 16S rRNA sequence; a well list (isolate_well.txt) including all detected ASV sequences and taxonomy. 

    # Claculate tables and figures, 2m
    time identify_isolate.R \
        --input result/ASV_table.txt \
        --genus result/genus.count \
        --taxonomy result/taxonomy_8.txt \
        --output result/isolate

- "Rarefaction curve in `result/isolate_rare_curve.pdf/png`"
- "Well list in `result/isolate_well.txt`"
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

    # (Optional) One plate heatmap for each library, delete '|head -n3' to plot all plates
    mkdir -p result/purity
    for l in `cat result/library.txt`; do
    # Prepare input file
    awk 'NR==FNR{a[$1]=$3}NR>FNR{print $1"\t"a[$1]}' result/isolate_well.txt seq/${l}.txt | \
        tail -n+2 | sed 's/\t$/\t0/' \
        > result/purity/${l}.txt
    # Format list into plate format
    format_list2plate.pl -i result/purity/${l}.txt \
        -o result/purity/${l}/
    # Batch plot top 1 plate, deleted '|head -n1' to plot all
    list=`ls result/purity/${l}/|cut -f 1 -d '.'|cut -f 2 -d 'P'|sort|uniq|head -n1`
    for plate in $(echo $list);do 
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
    graphlan_plot.sh -i $db -o result/graphlan

- result in: result/graphlan/graphlan.pdf

![image](http://210.75.224.110/github/Culturome/script/fig/graphlan.png)

**Fig 8. Cladogram showing the taxonomic distribution and occurrence frequency of cultivated bacteria.** The inner ring represents the dereplicated ASVs from cultivated root bacteria. Heat map in the outer ring represent the log2 transformed number of cultivated bacterial isolates belong to the corresponding ASV.

### 13. Project cleanup

When the project is complete, you can compress input files and delete temporary files to save space.

    # Stat disk usage
    du -sh *
    # Compress raw data and database
    gzip seq/*.fq
    # Delete temporary files
    rm -r temp
    # Packaging results
    zip result.zip -r result/

All the important results in `result` folder.


## ? TROUBLESHOOTING

### Wget file incomplete

wget download and rename problem, result only 8kb.
Try download and rename seperate.

    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_f1.fq.gz
    mv CRR127980_f1.fq.gz L1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa3/CRA002517/CRR127980/CRR127980_r2.fq.gz
    mv CRR127980_r2.fq.gz L1_2.fq.gz
    gunzip *.gz
    
## Change log

### 2020-08-12 v1.0

Finished single and multiple libraries analysis pipeline, include script, virtualbox and webserver

### 2022-03-21 v1.1

1. Make a conda package for personal computer run in Windows subsystems for Linux
2. Test on 20 GB data included 7 libraries, and record the times
3. vsearch replace usearch in cut primer step
4. Rscript replace stat of usearch
5. Test on WSL Ubuntu 20.04

### 2022-04-11 v2.0

1. Make a conda package for personal computer run in Mac 10.0+
2. Extensive testing, running successfully on over 50 computers
3. Add the Chinese pipeline 1pipeline.sh
4. Record the install and runing video, and provide the document in PPTX and DOCX
5. Open a QQ group for long-term supporint this project for backup files and discussion
6. Integrate single library and multi libraries pipeline into one
