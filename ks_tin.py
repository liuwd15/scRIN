#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
This program calculates KS, TIN and some other metrics for each transcript in BED file.

Requirements:
python >= 2.7
matplotlib >= 2.2.0
pysam >= 0.9
RseQC >= 2.6.4

Copyright (C) 2018 Wendao Liu
Copyright (C) 2016 original author Liguo Wang

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
-------------------------------------------------------------------------------------------------'''
import sys,os
import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from optparse import OptionParser
from qcmodule import getBamFiles
from time import strftime
from scipy.special import smirnov
from scipy.stats import rankdata


def printlog (message):
    '''
    message: string.
    Print message into stderr with time string appending to it.
    '''
    message="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + message
    print >>sys.stderr,message
    
def count_line(file):
    '''
    file: string, path of the file.
    Count the number of lines in the bedfile.
    '''
    count = 0  
    thefile=open(file,'rb')  
    while True:  
        buffer=thefile.read(1024*8192)  
        if not buffer:  
            break  
        count+=buffer.count('\n')  
    thefile.close()  
    return count 

def genomic_positions(refbed):
    '''
    refbed: string, reference BED file including all genes.
    exon: bool, determining whether to calculate metrics on exons rather than transcripts.
    Return genomic positions of each nucleotide in transcripts/exons.
    '''
    if refbed is None:
        print >>sys.stderr,"You must specify a bed file representing gene model\n"
        exit(0)
    
    for line in open(refbed,'r'):
        try:
            if line.startswith(('#','track','browser')):continue  
            #Parse fields from BED file. 
            fields = line.split()
            chrom = fields[0]
            chrom_start = int( fields[1] )
            name = fields[3]
            strand = fields[5]               
            exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
            exon_starts = map((lambda x: x + chrom_start ), exon_starts)
            exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
            exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
        except:
            print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
            continue
            
        #Return genomic positions of each nucleotide in transcripts.
        else:
            yield (name, chrom, exon_starts, exon_ends, strand)
        

def calculate_coverage(samfile, chrom, exon_starts, exon_ends):
    '''
    samfile: pysam.Samfile, input BAM file.
    chrom: string, index of chromosome.
    positions: list of integers, genomic positions of nucleotides.
    Return coverage on positions, a list of integers.
           coverage rate, a float number.
    '''
    #Calculate coverage of exons.
    if options.exon:
        exon_coverages = []
        covered_length = 0.0
        total_length = 0
        for (start, end) in zip(exon_starts, exon_ends):
            exon_coverage = np.array(samfile.count_coverage(chrom, start, end)).sum(axis=0)
            exon_coverages.append(exon_coverage)
            total_length += exon_coverage.size
            covered_length += sum(exon_coverage > 0)
        #Percentage of covered nucleotides.
        coverage_rate = covered_length / total_length
        return exon_coverages,coverage_rate
        
    #Calculate coverage of transcript.
    coverage = np.array([]) 
    for (start, end) in zip(exon_starts, exon_ends):
        exon_coverage = np.array(samfile.count_coverage(chrom, start, end)).sum(axis=0)
        coverage = np.append(coverage, exon_coverage)
    #Percentage of covered nucleotides.    
    coverage_rate = float(sum(coverage > 0)) /  coverage.size      
    return coverage,coverage_rate
    
def calculate_ks(coverage,read_length):
    '''
    coverage: array of integers, coverage on nucleotides.
    read_length: integer, the length of reads.
    Return KS and corresponding P value.
    '''
    #Calculate mean KS of exons.
    if options.exon:
        depth = np.concatenate(coverage).sum()
        mean_ks = 0.0
        for coverage_i in coverage:
            depth_i = coverage_i.sum()
            if depth_i == 0: continue
            length_i = coverage_i.size
            max_difference = 0.0 
            cumulative_coverage = 0.0
            for i in range(length_i): 
                cumulative_coverage += coverage_i[i]
                if abs(max_difference) < abs(cumulative_coverage/depth_i - float(i+1)/length_i):
                    max_difference = cumulative_coverage/depth_i - float(i+1)/length_i
            mean_ks += max_difference * depth_i / depth
        mean_p = smirnov(int(depth/read_length),abs(mean_ks))
        return mean_ks, mean_p
        
    #Calculate KS of transcript.
    length = coverage.size
    depth = np.sum(coverage)            #depth = read number * read length
    max_difference = 0.0 
    cumulative_coverage = 0.0
    for i in range(length):
        #Coverage here becomes empirical CDF of coverage.
        cumulative_coverage += coverage[i]
        if abs(max_difference) < abs(cumulative_coverage/depth - float(i+1)/length):
            max_difference = cumulative_coverage/depth - float(i+1)/length
    
    #Trapezoid uniform coverage.        
    #coverage_expect = []
    #if(read_length >= length):
    #    coverage_expect = [1/length] * int(length)
    #elif(2 * read_length > length):
    #    for k in range(int(length)):
    #        coverage_expect.append((min(k + 1, length - k, length - read_length + 1)) / (length - read_length + 1) / read_length)
    #else:
    #    for k in range(int(length)):
    #        coverage_expect.append((min(k + 1, length - k, read_length)) / (length - read_length + 1) / read_length)
    #d = 0
    #for x in range(int(length)):
    #    if abs(d) < abs(coverage[x]/depth - sum(coverage_expect[0:x])):
    #        d = coverage[x]/depth - sum(coverage_expect[0:x])
    
    ks = max_difference 
    p = smirnov(int(depth/read_length),abs(max_difference))
    return ks,p

def calculate_tin(coverage):
    '''
    coverage: array of integers, coverage on nucleotides.
    Return TIN value.
    '''
    #Calculate mean TIN of exons.
    if options.exon:
        length = np.concatenate(coverage).size
        depth = np.concatenate(coverage).sum()
        mean_tin = 0.0
        for coverage_i in coverage:            
            depth_i = coverage_i.sum()
            length_i = coverage_i.size
            entropy = 0.0
            for i in coverage_i: 
                if i > 0:
                    entropy += (float(i)/depth_i) * np.log(float(i)/depth_i)    
            mean_tin += 100*(np.exp(-entropy)) / length_i * depth_i / depth
        return mean_tin
    
    #Calculate TIN of transcript.
    length = coverage.size
    depth = coverage.sum()
    entropy = 0.0
    for i in coverage:
        if i > 0:
            entropy += (float(i)/depth) * np.log(float(i)/depth)    
    tin = 100*(np.exp(-entropy)) / length
    return tin

def mean_coverage(coverage):
    '''
    Calculate the mean coverage on a trancript.
    '''
    if options.exon:
        return np.mean(np.concatenate(coverage))
    
    return np.mean(coverage)

def coverage_size(coverage):
    '''
    Calculate the length of a transcript.
    '''
    if options.exon:
        return np.concatenate(coverage).size
    
    return coverage.size

def process_one_bamfile(file, processed_file):
    '''
    file: string, path of bamfile to be processed.
    processed_file: int, number of processed files.
    Write a metrics file in the same directory as bamfile.
    '''
    print >>sys.stderr, "\t" + file         
    printlog("Processing " + file)
    samfile = pysam.Samfile(file, "rb")
    #Details of each transcript/exon.
    OUT = open(file.replace('bam','') + 'metrics.xls','w')
    print >>OUT, "\t".join(["geneID","chrom", "length","average_coverage","coverage_rate","tin","ks","p"])
       
    processed_transcript = 0            #number of processed transcript/exon.
    read_length = samfile.head(1).next().infer_query_length()
    
    #Process each transcript/exon.
    for name, chrom, exon_starts, exon_ends, strand in genomic_positions(options.ref_bed):
        
        coverage,rate = calculate_coverage(samfile, chrom, exon_starts, exon_ends)
        
        #Omit transcripts/exons with average depth.
        if mean_coverage(coverage) < options.minimum_average_depth:
            ks_array[processed_file][processed_transcript] = 0
            tin_array[processed_file][processed_transcript] = 0
            print >>OUT, '\t'.join([str(i) for i in (name, chrom, coverage_size(coverage), mean_coverage(coverage), rate,  0, 0, 1)])
        
        #Process transcripts/exons above minimum average depth.            
        else:
            ks,p = calculate_ks(coverage,read_length)
            tin = calculate_tin(coverage)
            #This step make sure that 3' bias results in positive KS and 5' bias results in negative KS.
            if strand == '+':
                ks = -ks
            ks_array[processed_file][processed_transcript] = ks
            tin_array[processed_file][processed_transcript] = tin
            print >>OUT, '\t'.join([str(i) for i in (name, chrom, coverage_size(coverage), mean_coverage(coverage), rate,  tin, ks,  p)])
        processed_transcript += 1
        print >>sys.stderr, " %d/%d transcripts finished\r" % (processed_transcript, transcript_number),
    
    OUT.close()   
    samfile.close()

def process_bamfiles(files):
    '''
    files: list of string, bamfiles to be processed.
    Write a sample summary file in the current directory.
    '''
    #Summary file of all samples.
    SUM_SAMPLE = open('summary_sample.xls','w')
    print >>SUM_SAMPLE, "\t".join(['Bam_file','ks(mean)','ks(stdev)','tin(median)','tin(stdev)'])
    
    sample_names = []
    sample_median_tins = []
    sample_tin_stds = []
    sample_mean_kss = []

    #Process each bamfile.
    processed_file = 0
    for file in bamfiles:
        process_one_bamfile(file, processed_file)
        #Only sum transcripts with high expression for each sample.
        high_expression_transcripts = tin_array[processed_file] > 0
        print >>SUM_SAMPLE, "\t".join( [str(i) for i in (os.path.basename(file),\
                                        np.mean(ks_array[processed_file,high_expression_transcripts]),\
                                        np.std(ks_array[processed_file,high_expression_transcripts]),\
                                        np.median(tin_array[processed_file,high_expression_transcripts]),\
                                        np.std(tin_array[processed_file,high_expression_transcripts]))])
        sample_names.append(file.split('/')[-1][:-4])
        sample_median_tins.append(np.median(tin_array[processed_file,high_expression_transcripts]))
        sample_tin_stds.append(np.std(tin_array[processed_file,high_expression_transcripts]))
        sample_mean_kss.append(np.mean(ks_array[processed_file,high_expression_transcripts]))
        processed_file += 1
    
    SUM_SAMPLE.close()
    #Plot the median TIN of each sample.
    sample_number = len(sample_names)
    index = np.arange(sample_number)
    if sample_number > 50:
        plt.figure(figsize = [sample_number/6.5, sample_number/5])
    plt.bar(index, sample_median_tins, yerr=sample_tin_stds, xerr=[sample_mean_kss,[0]*sample_number])
    plt.ylabel('nonuniform coverage <-- median TIN --> uniform coverage')    
    plt.xticks(index, sample_names, rotation=90, fontsize=5)   
    plt.title('median TIN of samples')
    plt.savefig('sample_median_TIN.pdf')
    plt.clf()
    
def summary_transcript():
    '''
    Write a transcript summary file in the current directory.
    '''
    #Summary file of high expression transcripts.
    SUM_TRANSCRIPT = open('summary_transcript.xls','w')
    print >>SUM_TRANSCRIPT, "\t".join(["geneID","chrom",'ks(mean)','ks(stdev)','tin(median)','tin(stdev)'])
    processed_transcript = 0
    tins = []
    kss = [] 
    for line in open(options.ref_bed,'r'):
        if line.startswith(('#','track','browser')):continue  
        #Parse fields from BED file. 
        fields = line.split()
        chrom = fields[0]
        name = fields[3]
        #Only output transcript with high expression in more than 50% samples.
        if np.mean(tin_array[:,processed_transcript]>0) > 0.5:
            print >>SUM_TRANSCRIPT, "\t".join( [str(i) for i in (name, chrom,\
                                                np.mean(ks_array[:,processed_transcript]),\
                                                np.std(ks_array[:,processed_transcript]),\
                                                np.median(tin_array[:,processed_transcript]),\
                                                np.std(tin_array[:,processed_transcript]))])
            tins.append(np.median(tin_array[:,processed_transcript]))
            kss.append(np.mean(ks_array[:,processed_transcript]))
        processed_transcript += 1
        
    SUM_TRANSCRIPT.close()
    #Plot the mean KS and median TIN of transcripts.
    f = plt.figure()
    ax = f.add_subplot(111)
    plt.scatter(kss, tins, s=3)
    plt.ylabel('nonuniform coverage <-- median TIN --> uniform coverage')
    plt.xlabel("5' bias <-- mean KS --> 3' bias")
    plt.title('KS and TIN of transcripts')
    cor = round(np.corrcoef(kss, tins)[0,1],2)
    plt.text(0.1, 0.9,"r="+str(cor), ha='center', va='center', transform=ax.transAxes)
    plt.savefig('transcript_KS_and_TIN.pdf')
    plt.clf()
    
def output_rank():
    '''
    Write a file reporting TIN rank in all samples of transcripts. 
    '''
    rank_array = np.apply_along_axis(rankdata, 0, tin_array)
    TIN_RANK = open('TIN_rank.xls','w')
    print >>TIN_RANK, "\t".join(['gene'] + [i.split('/')[-1][:-4] for i in bamfiles])
    processed_transcript = 0
    for line in open(options.ref_bed,'r'):
        if line.startswith(('#','track','browser')):continue  
        #Parse fields from BED file. 
        fields = line.split()
        name = fields[3]
        if np.all(tin_array[:,processed_transcript]>0):
            print >>TIN_RANK, "\t".join([name] + [str(i) for i in rank_array[:,processed_transcript]])
        processed_transcript += 1
    TIN_RANK.close()
        

#Set options.
usage="%prog [options]" + '\n' + __doc__ + '\n'
parser = OptionParser(usage)
parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input BAM file(s). "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files (no spaces allowed). 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam files (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools. [required]')
parser.add_option("-r","--refbed",action="store",type="string",dest="ref_bed",help='Reference gene model in BED format. Must be strandard 12-column BED file. [required]')
parser.add_option("-d","--minDepth",action="store",type="int",dest="minimum_average_depth",default=1,help="Minimum average depth on each transcript/exon. default=%default")
parser.add_option("-e","--exon",action="store_true",dest="exon",help="Calculate adjusted KS and TIN on exons to reduce the effect of alternative splicing. This will mostly result in KS closer to 0 and bigger TIN.")
parser.add_option("-k","--rank",action="store_true",dest="rank",help="Output a xls file reporting the rank of transcript TIN across samples.")
(options,args)=parser.parse_args()

if not (options.input_files and options.ref_bed):
    parser.print_help()
    sys.exit(0)
    
if not os.path.exists(options.ref_bed):
    print >>sys.stderr, '\n\n' + options.ref_bed + " does NOT exists" + '\n'
    parser.print_help()
    sys.exit(0)

#Read BAM files.    
printlog("Get BAM file ...")
bamfiles = getBamFiles.get_bam_files(options.input_files)
transcript_number = count_line(options.ref_bed)
bamfile_number = len(bamfiles)
    
if bamfile_number <= 0:
    print >>sys.stderr, "No BAM file found, exit."
    sys.exit(0)
else:
    print >>sys.stderr, "Total %d BAM file(s):" % len(bamfiles)
        
#Matrix of ks and tin of all transcript in all samples.
ks_array = np.zeros((bamfile_number, transcript_number))
tin_array = np.zeros((bamfile_number, transcript_number))

process_bamfiles(bamfiles)
summary_transcript()
if options.rank:
    output_rank()