# Documentation
This program calculates KS, TIN and some other metrics for each transcript in BED file.

## Requirements:
You have to install or update some python packages before running this program.
python >= 2.7
matplotlib >= 2.2.0
pysam >= 0.9
RseQC >= 2.6.4
numpy any version
scipy any version

## Installation
This program is a python script and works in unix-like operating systems.
You can download it from [Github](https://github.com/liuwd15/Test/blob/master/ks_tin.py).
> wget https://github.com/liuwd15/Test/blob/master/ks_tin.py

Then you can run it with python.
> python ks_tin.py

Or you can make it executable and move it to your PATH.
> chmod +x ks_tin.py
> mv ks_tin.py ~/bin #"~/bin" can be replaced with other path in you PATH.
> ks_tin.py

## Usage
The required input files include:

* Sorted and indexed .bam file(s). You can use samtools to sort and index a .bam file.
> samtools sort -o example_sorted.bam example.bam
> samtools index example_sorted.bam

* Reference .bed file containing a list of gene models. Refseq transcripts are recommended.

The simplest usage is:
> ks_tin.py -r example_refseq.bed -i example_sorted.bam

To reduce the effect of alternative splicing in single cells, use **-e** option to calculate TIN on exons.
> ks_tin.py -r example_refseq.bed -i example_sorted.bam -e

The gene models with low expression are filted out. The default threshold is average read coverage (read length * mapped read number / gene model length) > 1, you can change it with **-d** option.
> ks_tin.py -r example_refseq.bed -i example_sorted.bam -d 2

It is also recommended that many .bam files should be processed together.
You can input comma-separated .bam files like this.
>ks_tin.py -r example_refseq.bed -i example1_sorted.bam,example2_sorted.bam,example3_sorted.bam

You can also create a text file containing the path of all .bam files like this.
> cat samples.txt
> ~/data/examples/example1/STAR_out/example1_sorted.bam
> ~/data/examples/example2/STAR_out/example2_sorted.bam
> ~/data/examples/example3/STAR_out/example3_sorted.bam

Then input the text file.
> ks_tin.py -r example_refseq.bed -i samples.txt

## Output
For each sample(.bam file), a .xls file containing some mRNA integrity metrics of each gene model will be created in the **same directory as .bam file**.
For all samples, a .xls file containing 2 overall mRNA integrity metrics (KS and TIN) and a .pdf file plotting overall mRNA integrity of each sample will be created in the **current directory**.
For all gene models, a .xls file containing 2 overall mRNA integrity metrics (KS and TIN) and a .pdf file plotting overall mRNA integrity of each gene model will be created in the **current directory**.
