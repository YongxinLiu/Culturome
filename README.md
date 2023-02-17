[TOC]

# High-throughput cultivation and identification of bacteria from the microbiota (Culturome)

This is the computational part of the detailed protocol for high-throughput bacterial isolation and cultivation. Now is published in Nature Protocols.

> Jingying Zhang, [Yong-Xin Liu](http://bailab.genetics.ac.cn/YongxinLiuEn.html), Xiaoxuan Guo, Yuan Qin, Ruben Garrido-Oter, Paul Schulze-Lefert & Yang Bai. (2021). High-throughput cultivation and identification of bacteria from the plant root microbiota. ***Nature Protocols***, 16(2): 988-1012: https://doi.org/10.1038/s41596-020-00444-7 (Highly Cited Paper) <span class="__dimensions_badge_embed__" data-doi="10.1038/s41596-020-00444-7" data-hide-zero-citations="true" data-style="small_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>


The overview of the workflow is in the following figures. 


![image](http://bailab.genetics.ac.cn/github/Culturome/script/fig/fig1.jpg)

**Figure 1 | Overview of the high-throughput bacterial cultivation and identification system.**

![image](http://bailab.genetics.ac.cn/github/Culturome/script/fig/fig4.jpg)

**Figure 2 | Workflow of bioinformatic analysis to identify cultivated bacteria**. a, Summary of the steps from obtaining raw data to generating clean sequences of the amplified 16S rRNA gene region. b, Workflow to identify amplicon sequence variants (ASVs) and taxonomy. c, Sequence and taxonomy information for cultivated bacteria in 96-well cell culture plates.

![image](http://bailab.genetics.ac.cn/github/Culturome/script/fig/fig6.jpg)

**Fig. 3 | Anticipated results for bacterial cultivation.**
a. Bar plot showing the distribution of read counts in each plate. b. Read counts for the negative (nuclease-free water as the PCR template) and positive controls (*E. coli* DNA as the PCR template). c. Rarefaction curve of ASVs based on the number of wells containing bacteria. d. Distribution of wells regarding the purity based on evaluation of ASVs of the bacterial 16S rRNA gene. e. An example of the purity of wells containing cultivated bacteria in a 96-well cell culture plate. f. Cladogram showing the taxonomic distribution and occurrence frequency of cultivated bacteria.

## Equipment

Personal computer running with 64-bit Linux Ubuntu 16.04+ or CentOS 7.5+, at least dual-core 2.4 GHz, 4 GB of RAM (8 GB preferred) and 30 GB space. The analysis pipeline is mainly based on QIIME, VSEARCH, GraPhlAn and home-made scripts (Culturome pipeline). All the home-made scripts and manual is deposit in https://github.com/YongxinLiu/Culturome. The dependency databases is RDP. To simplify the installation procedure for all computational software, we provide a conda pack for download. Additionally, a graphic user interface Webserver have be construct for user to analysis culturomics data online http://bailab.genetics.ac.cn/culturome/, and fully installed VirtualBox image http://bailab.genetics.ac.cn/culturome/download/QIIME-1.9.1-amd64.zip (6 GB).

### Computer requirements

- Operational system requirements: 64-bit Linux, such as Ubuntu 16.04+ or CentOS 7.5+, Windows Subsystem for Linux Ubuntu 20.04 LTS
- Hardware requirements: > 4 GB+ RAM and 10 GB free disk space.
- Hardware recommended: > 8 GB+ RAM, 30GB free disk

### Software

- Conda 4.8.3+, software management system (https://repo.anaconda.com/miniconda/).
- R 4.0+ (https://cran.r-project.org). The analysis was tested using R 4.1. Packages include ggplot2, pheatmap, dplyr, stringr should be installed. 
- QIIME version 1.9.1 http://qiime.org/ for split barcodes
- VSEARCH v.2.7.1 or later https://github.com/torognes/vsearch
- GraPhlAn 1.1.3-1, for visualizing taxonomic tree
- Culturome v2.0 analysis pipeline, related scripts, and example results are available at https://github.com/YongxinLiu/Culturome

### Input data files

- Sequencing Data：The raw sequence data used in this paper have been deposited in the Genome Sequence Archive (GSA) in BIG Data Center, under accession numbers CRA002517. It can be download automatic in running demo pipeline. 
- RDP train set 16 as default taxonomy database is include in Culturome v2.0. The more taxonomy databases are in http://www.drive5.com/sintax.

## Equipment setup

For Linux system, conda (https://docs.conda.io/en/latest/miniconda.html) is importantance to setup the pipeline. We tested the pipeline on Ubuntu 20.04.

To install Conda, download the latest version from https://docs.conda.io/en/latest/miniconda.html, install by bash, and configure the Bioconda58 channel for easy-install software in biology: 

For Linux users, include Windows subsystem for Linux

    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    ~/miniconda3/condabin/conda init

For Windows user, install R, RStuido and Windows Subsystem for Linux (WSL) are requirement, and the full install steps in pipeline.md (English version) and 1pipeline.sh 附录1 (中文版)

For Mac users, install xcode, homebrew, and coreutils are requirement,  and the full install steps in pipeline.md (English version) and 1pipeline.sh 附录2 (中文版)

For the backup file, and discuss. Please add QQ group chat (QQ群：815816627).

### Method 1. Quick download (Recommended)

For Linux users, include Windows subsystem for Linux

    # Install culturome conda environment
    n=culturome
    wget -c http://bailab.genetics.ac.cn/db/conda/${n}.tar.gz
    mkdir -p ~/miniconda3/envs/${n}
    tar -xvzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    conda activate $n
    conda unpack
    rm culturome.tar.gz

    # Install culturome pipeline
    # For windows usea change to c:
    wd=/mnt/c/microbiome 
    mkdir -p $wd && cd $wd
    wget -c http://bailab.genetics.ac.cn/db/Culturome.tar.gz
    tar xvzf Culturome.tar.gz
    chmod +x Culturome/script/*
    echo export PATH=`pwd`/Culturome/script:\$PATH >>  ~/.bashrc
    # Close the terminal, and reopen    

For Mac users

    n=culturome
    wget -c http://bailab.genetics.ac.cn/db/conda/${n}_mac.tar.gz
    mv ${n}_mac.tar.gz ${n}.tar.gz
    mkdir -p ~/miniconda3/envs/${n}
    tar -xvzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    conda activate $n
    conda unpack
    
    # Install culturome pipeline
    wd=~/Documents/microbiome
    mkdir $wd && cd $wd
    wget -c http://bailab.genetics.ac.cn/db/Culturome.tar.gz
    tar xvzf Culturome.tar.gz
    cp -f Culturome/mac/* Culturome/script/
    chmod +x Culturome/script/*
    echo export PATH=`pwd`/Culturome/script:\$PATH >> ~/.zshrc
    source ~/.zshrc

### Method 2. Install pre-requisite by conda

According the tip: type `yes` to bash Miniconda2-latest-Linux-x86_64.sh, then Press `ENTER` to confirm the location. After install finished, close and re-open your current terminal to effect the changes.

    conda config --add channels conda-forge
    conda config --add channels bioconda

(Optinal) For Chinese user only, add tsinghua mirror to accelerate download.

    site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
    conda config --add channels ${site}/cloud/conda-forge/
    conda config --add channels ${site}/cloud/bioconda/

Create a new environment culturome
    
	cconda create -n culturome python=2.7 -c bioconda -y 
	conda activate culturome
	conda install qiime=1.9.1 -c bioconda
    conda install graphlan=1.1.3 -c bioconda
    conda install r-base r-dplyr r-optparse r-stringr r-ggplot2 r-pheatmap r-reshape2 -c bioconda
    
More install methods or troubleshooting in http://qiime.org/install/install.html

If you install with problem, a QIIME VirtualBox image can be downloaded from http://qiime.org/install/virtual_box.html, and run in Windows/Mac/Linux system following the instructions.

To install Cultrome pipeline, download them from Github, unzip the zip file and add the directory to your PATH environment variable: 

    # Method1. Download from github
    wget -c https://github.com/YongxinLiu/Culturome/archive/master.zip
    unzip master.zip
    mv Culturome-master Culturome
    
    # Method2. Download from backup site
    wget -c http://bailab.genetics.ac.cn/db/Culturome.tar.gz
    tar xvzf Culturome.tar.gz
    
    # Add environment variables permanently 
    chmod +x Culturome/script/*
    echo export PATH=`pwd`/Culturome/script:\$PATH >>  ~/.bashrc
    # Then close the terminal, and restart to effect

### Method 3. VirtualBox & Webserver

To simplify the installation procedure for all computational software, we provide a fully installed VirtualBox image http://bailab.genetics.ac.cn/culturome/download/QIIME-1.9.1-amd64.zip, and a webserver https://github.com/YongxinLiu/Culturome.

## Analysis protocols

Using Rstudio open the 1pipeline.sh, then start the analysis step by step.

## Contact

Prof. Yong-Xin Liu

Yong-Xin Liu's lab, Agricultural Genomics Institute at Shenzhen, Chinese Academy of Agricultural Sciences

No 97, Buxin Road, Dapeng District, Shenzhen 518120, Guangdong, China

E-mail: liuyongxin@caas.cn

Wechat: meta-genomics

Cite: Jingying Zhang, [Yong-Xin Liu](https://agis.caas.cn/en/research/principalinvestigator/4510ae674c0844289e50ae538d4bb6c3.htm), Xiaoxuan Guo, Yuan Qin, Ruben Garrido-Oter, Paul Schulze-Lefert & Yang Bai. (2021). High-throughput cultivation and identification of bacteria from the plant root microbiota. ***Nature Protocols***,  https://doi.org/10.1038/s41596-020-00444-7 <span class="__dimensions_badge_embed__" data-doi="10.1038/s41596-020-00444-7" data-hide-zero-citations="true" data-style="small_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
