﻿# scRIN Documentation

---


## Introduction
This python script is used for measuring the mRNA integrity with single-cell sequencing data. The analysis is conducted on 3 levels, sample/cell, gene/transcript and exon.    
mRNA integrity is measured by 2 criteria, KS and TIN. KS measures the read coverage bias and TIN measures the read coverage uniformity on each gene model.

## Requirements
You have to install or update some python packages before running this program.   

python >= 2.7   
matplotlib >= 2.2.0   
pysam >= 0.9   
RSeQC >= 2.6.4   
numpy  

## Installation
This program is a python script and works in unix-like operating systems.
You can download it from [Github](https://raw.githubusercontent.com/liuwd15/scRIN/master/scRIN.py).

    wget https://raw.githubusercontent.com/liuwd15/scRIN/master/scRIN.py

Then run it with python.

    python scRIN.py

Or you can make it executable and move it to your PATH.

    chmod +x scRIN.py  
    mv scRIN.py \~/bin #"\~/bin" can be replaced with other path in you PATH.  
    scRIN.py

## Usage
### Required input

1. Sorted and indexed **.bam** file(s). Samtools can be used to sort and index a .bam file.

        samtools sort -o example_sorted.bam example.bam  
        samtools index example_sorted.bam

2. Reference [12-column](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) **.bed** file containing a list of gene models. Representative .bed file containing RefSeq transcripts of hg19 and mm10 are available [here](https://github.com/liuwd15/scRIN/tree/master/bed). Only the longest transcript for the gene with multiple transcripts is included to avoid redundancy.

The simplest usage is:

    scRIN.py -r example_refseq.bed -i example_sorted.bam

It is also recommended that many .bam files should be processed together.
You can input comma-separated .bam files like this.

    scRIN.py -r example_refseq.bed -i example1_sorted.bam,example2_sorted.bam,example3_sorted.bam

You can also create a text file containing the path of all .bam files like this.

    cat samples.txt  
    ~/data/examples/example1/example1_sorted.bam  
    ~/data/examples/example2/example2_sorted.bam  
    ~/data/examples/example3/example3_sorted.bam

Then input the text file.

    scRIN.py -r example_refseq.bed -i samples.txt

### Other options

To perform additional analysis on inter-exon, use **-e** option. This will result in more output files (see Output).

    scRIN.py -r example_refseq.bed -i example_sorted.bam -e

The transcripts with low expression are filted out. The default threshold is average coverage (read length * mapped read number / gene model length) > 0.5, you can change it with **-d** option.

    scRIN.py -r example_refseq.bed -i example_sorted.bam -d 1
    
By default, only transcript/exon expressed in more than 3 samples will be summarized in **summary_transcript.xls** file (see Output), you can change this threshold with **-s** option. Noticing that only when the number of input files is more than this threshold will transcript/exon be summarized.

    scRIN.py -r example_refseq.bed -i example_sorted.bam -s 5

If you want to get the rank of TIN of each trancript across samples, use **-k** option. This will create .xls files containing the rank of TIN of transcripts expresses in all samples.

    scRIN.py -r example_refseq.bed -i example_sorted.bam -k

## Output
### Sample directory
For each sample(.bam file), following files will be created in the same directory as .bam file.

1. A **.metrics.xls** file containing some mRNA integrity metrics of each transcript. Metrics include:

* `average_coverage`: read length * mapped read number / gene model length
* `coverage_rate`: length of read mapped region / total length
* `exon_number`: number of exons
* `base_level_KS`: Measurment of read coverage bias on the whole transcript. Ranging from -1 to 1. Value close to -1 suggests 5' bias of read coverage, and value close to 1 suggests 3' bias.
* `inter_exon_KS`: Measurment of read coverage bias between exons. Ranging from -1 to 1. Value close to -1 suggests read counts on exons near 5' end are generally bigger than those near 3' end, and value close to 1 suggests the opposite.
* `intra_exon_KS`: Measurment of read coverage bias within exons. Ranging from -1 to 1. Value close to -1 suggests read counts within exons generally have 5' bias, and value close to 1 suggests the opposite.
* `base_level_TIN`: Measurement of read coverage uniformity on the whole transcript. Ranging from 0 to 100, and high value suggests strong uniformity.
* `inter_exon_TIN`: Measurement of read coverage uniformity between exons. Ranging from 0 to 100, and high value suggests strong uniformity.
* `intra_exon_TIN`: Measurment of read coverage uniformity within exons. Ranging from 0 to 100, and high value suggests strong uniformity.

2. A **.KS_TIN.pdf** file containing 3(*4*) scatter plots.

* `base_level_TIN` vs `base_level_KS` for each transcript.
* `inter_exon_TIN` vs `inter_exon_KS` for each transcript.
* `intra_exon_TIN` vs `intra_exon_KS` for each transcript.
* *[option -e enabled]* `exon_TIN` vs `exon_KS` for each exon.

3. *[option -e enabled]* A **.exon.xls** file containing metrics of each exon. Metrics include:
* `exon_KS`: Measurment of read coverage bias on the exon.
* `exon_TIN`: Measurement of read coverage uniformity on the exons.

### Current directory
Following files will be created in the current directory.

1. A **summary_sample.xls** file containing the summary metrics of each sample. Metrics include:

* `expressed_transcript`: Number of expressed transcripts.
* `base_level_KS(mean)`: Mean of transcript KS.
* `inter_exon_KS(mean)`: Mean of inter-exon KS.
* `intra_exon_KS(mean)`: Mean of intra-exon KS.
* `base_level_TIN(mean)`: mean of transcript TIN.
* `inter_exon_TIN(mean)`: mean of inter-exon TIN.
* `intra_exon_TIN(mean)`: mean of intra-exon TIN.
* `base_level_KS(std)`: Standard deviation of transcript KS.
* `inter_exon_KS(std)`: Standard deviation of inter-exon KS.
* `intra_exon_KS(std)`: Standard deviation of intra-exon KS.
* `base_level_TIN(std)`: Standard deviation of transcript TIN.
* `inter_exon_TIN(std)`: Standard deviation of inter-exon TIN.
* `intra_exon_TIN(std)`: Standard deviation of intra-exon TIN.

2. A **summary_sample.pdf** file containing 4 barplots.
* `base_level_TIN(mean)` with `base_level_TIN(std)` as error bar for each sample.
* `inter_exon_TIN(mean)` with `inter_exon_TIN(std)` as error bar for each sample.
* `intra_exon_TIN(mean)` with `intra_exon_TIN(std)` as error bar for each sample.
* `expressed_transcript`.

3. A **summary_transcript.xls** file containing the summary metrics of transcript expressed in all samples. Metrics are similar as those in 1 but are calculated across samples.

4. A **summary_transcript.pdf** file containing 3 scatter plots.

* `base_level_TIN(mean)` vs `base_level_KS(mean)`
* `inter_exon_TIN(mean)` vs `inter_exon_KS(mean)`
* `intra_exon_TIN(mean)` vs `intra_exon_KS(mean)`

5. *[option -e enabled]* A **summary_exon.xls** file containing the summary metrics of exon expressed in all samples. Metrics include:

* `exon_KS(mean)`: Mean of exon KS across samples.
* `exon_TIN(mean)`: mean of exon TIN across samples.
* `exon_KS(std)`: Standard deviation of exon KS across samples.
* `exon_TIN(std)`: Standard deviation of exon TIN across samples.

6. *[option -k enabled]* Three files **base_level_TIN_rank.xls**, **inter_exon_TIN_rank.xls** and **intra_exon_TIN_rank.xls** containing the ranks of TIN of each transcript across all samples.