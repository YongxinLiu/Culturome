wd=/mnt/bai/yongxin/github/Culturome
PATH=${wd}/script:$PATH
mkdir -p seq temp result

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
    
# 2. Validate the mapping file, 1s
for l in `cat result/library.txt`; do
validate_mapping_file.py -m seq/${l}.txt \
    -o temp/
done

# 3. Merge pair-end reads, 5s
for l in `cat result/library.txt`; do
time vsearch -fastq_mergepairs seq/${l}_1.fq \
    -reverse seq/${l}_2.fq \
    -fastqout temp/${l}.fq
done

# 4. Demultiplexing

# extract barcodes, 1m37s, note barcode legnth
rm -rf temp/qc.fa
for l in `cat result/library.txt`; do
time extract_barcodes.py \
    -f temp/${l}.fq -m seq/${l}.txt \
	-c barcode_paired_stitched \
	--bc1_len 10 --bc2_len 6 \
	-a --rev_comp_bc2 \
	-o temp/${l}

# split library, 2m8s
time split_libraries_fastq.py \
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

# 5. Visualize counts of samples

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

# 6. Cut primers
# length of primers used, 10s
time usearch -fastx_truncate temp/qc.fa \
	-stripleft 19 -stripright 18 \
	-fastaout temp/filtered.fa

# 7. Pick representitve sequences
# Calculate frequency of non-reduncancy reads, 4s
time vsearch \
    --derep_fulllength temp/filtered.fa \
	--relabel Uni --minuniquesize 8 --sizeout \
	--output temp/uniques.fa 

# Denoise by unoise3, 2s
time usearch -unoise3 temp/uniques.fa \
    -zotus temp/Zotus.fa

# Rename to ASV
awk 'BEGIN {n=1}; />/ {print ">ASV_" n; n++} !/>/ {print}' temp/Zotus.fa \
    > temp/otus.fa
cp temp/otus.fa result/otu.fa

# 8. (Optional) Remove host and unspecific amplification

# 9. Construct ASV table
# 96 threads, 99.07% matched, 38s, 26m
time vsearch --usearch_global temp/filtered.fa \
    --db result/otu.fa \
    --otutabout temp/otutab.txt \
    --id 0.97

# Stat, e.g. 2678 samples, 378 OTUs
usearch -otutab_stats temp/otutab.txt \
    -output temp/otutab.stat
#cat temp/otutab.stat

# 10. Calculate 100%(1) false discovery reads, mean cut all negative control
negative_threshold.R \
    --input temp/otutab.txt \
    --metadata result/metadata.txt \
    --threshold 1 \
    --negative A12 \
    --positive B12 \
    --output result/fdr.txt

# Filter flase discovery well in feature table
usearch -otutab_trim temp/otutab.txt \
    -output result/otutab.txt \
    -min_sample_size `cat result/fdr.txt` 

# 11. Taxonomic classification
# RDP 4s, 10 threads
time usearch -sintax result/otu.fa \
    -db ${wd}/db/rdp_16s_v16_sp.fa \
	-tabbedout temp/otu.fa.tax \
	-sintax_cutoff 0.6 -strand both
tax_sum.sh -i temp/otu.fa.tax

# Prepare graphlan files
graphlan_prepare_order.R \
    --input result/otutab.txt \
    --output result/graphlan/ \
    --taxonomy result/taxonomy_8.txt \
    --abundance 0 \
    --number 1000
# Plot graphlan
graphlan_plot.sh -i ${wd} \
    -o result/graphlan

# 12. Identify non-redundancy isolates
# Input ASV table and taxonomy, and output tables and figure
time identify_isolate.R \
    --input result/otutab.txt \
    --taxonomy result/taxonomy_8.txt \
    --output result/isolate

# 13. Each well purity in plate
# Fetch purity in result table
mkdir -p result/purity
for l in `cat result/library.txt`; do
awk 'NR==FNR{a[$1]=$3}NR>FNR{print $1"\t"a[$1]}' result/isolate_sample.txt seq/${l}.txt | tail -n+2 | sed 's/\t$/\t0/' \
    > result/purity/${l}.txt
# Format list into plate format
format_list2plate.pl -i result/purity/${l}.txt \
    -o result/purity/${l}/

# Batch plot each plate
list=`ls result/purity/${l}/|cut -f 1 -d '.'|cut -f 2 -d 'P'|sort|uniq`
time for plate in $list;do 
    plot_pheatmap.sh -i result/purity/${l}/${l}P${plate}.plate \
        -o result/purity/${l}/${l}P${plate}.plate.heat.png
done
done