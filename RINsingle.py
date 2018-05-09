#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
This program calculates KS, TIN and some other metrics for each transcript in BED file.

Requirements:
python >= 2.7
matplotlib >= 2.2.0
pysam >= 0.9
RseQC >= 2.6.4

Copyright (C) 2018 Wendao Liu
Copyright (C) 2016 Liguo Wang, original author of some functions 

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
#from scipy.special import smirnov
from scipy.stats import rankdata


def printlog(message):
    '''
    Print message into stderr with time string appending to it.
    message: string.    
    '''
    message='@ ' + strftime('%Y-%m-%d %H:%M:%S') + ': ' + message
    print >>sys.stderr,message
    
def count_number(refbed):
    '''
    Count the number of transcripts and exons in the bedfile.
    file: string, path of the file. 
    Return transcript_number and exon_number, 2 integers.
    '''
    transcript_number = 0
    exon_number = 0
    if refbed is None:
        print >>sys.stderr,'You must specify a bed file representing gene model\n'
        exit(0)
    
    for line in open(refbed,'r'):
        try:
            if line.startswith(('#','track','browser')):continue  
            #Parse fields from BED file.   
            exon_number += int(line.split()[9])
            transcript_number += 1
        except:
            print >>sys.stderr, '[NOTE:input bed must be 12-column] skipped this line: ' + line,
            continue 
    return transcript_number, exon_number

def genomic_positions(refbed):
    '''
    Locate genomic positions of each nucleotide in exons of each transcript.
    refbed: string, reference BED file including all transcripts.
    Return the name, chromosome, number of exons, start and end position of exons and strand.
    '''
    if refbed is None:
        print >>sys.stderr,'You must specify a bed file representing gene model\n'
        exit(0)
    
    for line in open(refbed,'r'):
        try:
            if line.startswith(('#','track','browser')):continue  
            #Parse fields from BED file. 
            fields = line.split()
            chrom = fields[0]
            chrom_start = int(fields[1])
            name = fields[3]
            strand = fields[5]     
            exon_number = int(fields[9])
            exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
            exon_starts = map((lambda x: x + chrom_start ), exon_starts)
            exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
            exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
        except:
            print >>sys.stderr,'[NOTE:input bed must be 12-column] skipped this line: ' + line,
            continue
            
        #Return genomic positions of each nucleotide in transcripts.
        else:
            yield (name, chrom, exon_number, exon_starts, exon_ends, strand)
        

def calculate_coverage(samfile, chrom, exon_starts, exon_ends):
    '''
    Calculate the coverages on exons.
    samfile: pysam.Samfile, input BAM file.
    chrom: string, index of chromosome.
    exon_starts: list of integers, genomic positions of the start of the exons.
    exon_starts: list of integers, genomic positions of the end of the exons.
    Return exon_coverages, a list of array of integers.
           total_length, a float number.
           average_coverage, a float number.
           coverage_rate, a float number.           
    '''
    exon_coverages = []
    covered_length = 0.0
    total_length = 0           #Length of the transcript.
    for (start, end) in zip(exon_starts, exon_ends):
        exon_coverage = np.array(samfile.count_coverage(chrom, start, end)).sum(axis=0)
        exon_coverages.append(exon_coverage)
        total_length += exon_coverage.size
        covered_length += sum(exon_coverage > 0)   
    coverage_rate = covered_length / total_length          #Percentage of covered nucleotides.
    average_coverage = np.mean(np.concatenate(exon_coverages))       #Average coverage on trasncript.
    
    return exon_coverages, total_length, average_coverage, coverage_rate
        

def ks_value(coverage, strand):
    '''
    Calculate the KS value of the coverage.
    coverage: array of integers, coverage on one exon of transcript.   
    Return max_difference, a float number.
    '''    
    depth = coverage.sum()
    if depth == 0:
        return 0
    length = coverage.size
    max_difference = 0.0 
    cumulative_coverage = 0.0
    for i in range(length):
        #Coverage here becomes empirical CDF of coverage.
        cumulative_coverage += coverage[i]
        if abs(max_difference) < abs(cumulative_coverage/depth - float(i+1)/length):
            max_difference = cumulative_coverage/depth - float(i+1)/length
    
    if strand == '-':
        return max_difference 
    return -max_difference

def tin_value(coverage):
    '''
    Calculate the TIN value of the coverage.
    coverage: array of integers, coverage on one exon of transcript.   
    Return TIN, a float number.
    '''  
    depth = coverage.sum()
    #if depth == 0:
    #    return 100
    length = coverage.size
    entropy = 0.0
    for i in coverage: 
        if i > 0:
            entropy += (float(i)/depth) * np.log(float(i)/depth) 
            
    return 100*(np.exp(-entropy)) / length

def calculate_ks_tin(exon_coverages, name, strand, EXON, processed_file, processed_exon):
    '''
    Calculate several KS and TIN values of the transcript. 
    exon_coverages: list of array of integers, coverage on exons of one transcript.
    name: string, the name of the transcript.
    strand: string, either '+' or '-', the strand of the transcript.
    EXON: file, optional file to output intra-exon KS and TIN.
    processed_file: int, number of processed files.
    processed_exon: int, number of processes exons.
    Return two arrays: base-level KS, exon-level KS, exon-average KS and same for TIN.
    '''
    exon_kss = []
    exon_tins = []
    exon_average_coverages = []
    depth = np.concatenate(exon_coverages).sum()
    exon_average_ks = 0.0
    exon_average_tin = 0.0
    
    for exon_coverage in exon_coverages:
        exon_average_coverages.append(exon_coverage.mean())
        #This step make sure that 3' bias results in positive KS and 5' bias results in negative KS.     
        exon_ks = ks_value(exon_coverage, strand)
        exon_kss.append(exon_ks)
        exon_average_ks += exon_ks * exon_coverage.sum() / depth
        
        exon_tin = tin_value(exon_coverage)
        exon_tins.append(exon_tin)
        exon_average_tin += exon_tin * exon_coverage.sum() / depth
      
    #KS and TIN of full length transcript.    
    base_level_ks = ks_value(np.concatenate(exon_coverages), strand)   
    base_level_tin = tin_value(np.concatenate(exon_coverages))
    #KS and TIN calculated by averaging coverage on each exon.
    exon_level_ks = ks_value(np.array(exon_average_coverages), strand)
    exon_level_tin = tin_value(np.array(exon_average_coverages))  
    
    if options.exon:
        exon_number = len(exon_kss)
        exon_array[processed_file, processed_exon:(processed_exon+exon_number), 0] = exon_kss
        exon_array[processed_file, processed_exon:(processed_exon+exon_number), 1] = exon_tins
        for i in range(exon_number):
            print >> EXON, '\t'.join([name+'.'+str(i+1), str(exon_average_coverages[i]),\
                                      str(exon_kss[i]), str(exon_tins[i])])
               
    return np.array([base_level_ks, exon_level_ks, exon_average_ks]),\
           np.array([base_level_tin, exon_level_tin, exon_average_tin])


def process_one_bamfile(file, processed_file):
    '''
    Prcess one bamfile.
    file: string, path of bamfile to be processed.
    processed_file: int, number of processed files.
    Output a metrics text file and a plot pdf file in the same directory as bamfile.
    '''
    print >>sys.stderr, '\t' + file         
    printlog('Processing ' + file)
    samfile = pysam.Samfile(file, 'rb')
    #Details of each transcript.
    OUT = open(file.replace('bam','') + 'metrics.xls','w')
    print >>OUT, '\t'.join(['gene','chrom','length','exon_number','average_coverage','coverage_rate',\
                            'base_level_KS','exon_level_KS','exon_average_KS',\
                            'base_level_TIN','exon_level_TIN','exon_average_TIN'])
    EXON = ''
    if options.exon:
        EXON = open(file.replace('bam','') + 'exon_metrics.xls','w')
        print >>EXON, '\t'.join(['exon','average_coverage','exon_KS','exon_TIN'])
    
    processed_transcript = 0            #Number of processed transcript.
    processed_exon = 0
    #read_length = samfile.head(1).next().infer_query_length()
       
    #Process each transcript/exon.
    for name, chrom, exon_number, exon_starts, exon_ends, strand in genomic_positions(options.ref_bed):        
        coverages, length, average, rate = calculate_coverage(samfile, chrom, exon_starts, exon_ends)
        
        #Omit unexpressed transcripts.
        if average == 0:
            processed_transcript += 1
            processed_exon += exon_number
            print >>OUT, '\t'.join([str(i) for i in (name, chrom, length, exon_number, average, rate,\
                                    0, 0, 0, 0, 0, 0)])
            continue
        #Process expressed transcripts.
        else:
            three_ks,three_tin = calculate_ks_tin(coverages, name, strand, EXON, processed_file, processed_exon)        
            ks_array[processed_file,processed_transcript] = three_ks
            tin_array[processed_file,processed_transcript] = three_tin
            if average < options.minimum_average_depth:
                expression_level[processed_file,processed_transcript] = 1  
                print >>OUT, '\t'.join([str(i) for i in (name, chrom, length, exon_number, average, rate,\
                                        0, 0, 0, 0, 0, 0)])                             
            else:
                expression_level[processed_file,processed_transcript] = 2
                print >>OUT, '\t'.join([str(i) for i in (name, chrom, length, exon_number, average, rate,\
                                        three_ks[0], three_ks[1], three_ks[2],\
                                        three_tin[0], three_tin[1], three_tin[2])])
            
        processed_transcript += 1
        processed_exon += exon_number
        print >>sys.stderr, ' %d/%d transcripts finished\r' % (processed_transcript, transcript_number),
        
    OUT.close()   
    samfile.close()
    index = expression_level[processed_file] == 2       
    
    #Plot the KSs and TINs of transcripts.
    if options.exon:
        EXON.close()
        exon_index = exon_array[processed_file,:,0] != 0
        plt.figure(figsize=[20,16])
        plt.subplot(221)
        plt.scatter(ks_array[processed_file,index,0], tin_array[processed_file,index,0], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Base-level KS and TIN')
        plt.subplot(222)
        plt.scatter(exon_array[processed_file,exon_index,0], exon_array[processed_file,exon_index,1], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Intra-exon KS and TIN')
        plt.subplot(223)
        plt.scatter(ks_array[processed_file,index,1], tin_array[processed_file,index,1], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Exon-level KS and TIN')
        plt.subplot(224)
        plt.scatter(ks_array[processed_file,index,1], tin_array[processed_file,index,1], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Exon-average KS and TIN')
        plt.savefig(file.replace('bam','') + 'KS_TIN.pdf')
        plt.clf()
    else:
        plt.figure(figsize=[20,16])
        plt.subplot(221)
        plt.scatter(ks_array[processed_file,index,0], tin_array[processed_file,index,0], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Base-level KS and TIN')
        plt.subplot(223)
        plt.scatter(ks_array[processed_file,index,1], tin_array[processed_file,index,1], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Exon-level KS and TIN')
        plt.savefig(file.replace('bam','') + 'KS_TIN.pdf')
        plt.subplot(224)
        plt.scatter(ks_array[processed_file,index,1], tin_array[processed_file,index,1], s=3)
        plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
        plt.xlabel("5' bias <-- KS --> 3' bias")
        plt.title('Exon-average KS and TIN')
        plt.savefig(file.replace('bam','') + 'KS_TIN.pdf')
        plt.clf()

def process_bamfiles(files):
    '''
    Process all bamfiles.
    files: list of string, bamfiles to be processed.
    Output a sample summary file and a plot pdf file in the current directory.
    '''
    #Summary file of all samples.
    SUM_SAMPLE = open('summary_sample.xls','w')
    print >>SUM_SAMPLE, '\t'.join(['Bam_file','expressed_transcript',\
                                   'base_level_KS(mean)','exon_level_KS(mean)','exon_average_KS(mean)',\
                                   'base_level_TIN(median)','exon_level_TIN(median)','exon_average_TIN(median)',\
                                   'base_level_KS(sd)','exon_level_KS(sd)','exon_average_KS(sd)',\
                                   'base_level_TIN(sd)','exon_level_TIN(sd)','exon_average_TIN(sd)'])
    sample_names = []
    sample_median_tins = []
    sample_tin_sds = []
    expressed_transcript = []

    #Process each bamfile.
    processed_file = 0
    for file in bamfiles:
        process_one_bamfile(file, processed_file)
        #Only sum transcripts with high expression for each sample.
        index = expression_level[processed_file] == 2
        expressed = sum(expression_level[processed_file] > 0)
        expressed_transcript.append(expressed)
        print >>SUM_SAMPLE, '\t'.join( [str(i) for i in (os.path.basename(file),expressed,\
                                        np.mean(ks_array[processed_file,index,0]),\
                                        np.mean(ks_array[processed_file,index,1]),\
                                        np.mean(ks_array[processed_file,index,2]),\
                                        np.median(tin_array[processed_file,index,0]),\
                                        np.median(tin_array[processed_file,index,1]),\
                                        np.median(tin_array[processed_file,index,2]),\
                                        np.std(ks_array[processed_file,index,0]),\
                                        np.std(ks_array[processed_file,index,1]),\
                                        np.std(ks_array[processed_file,index,2]),\
                                        np.std(tin_array[processed_file,index,0]),\
                                        np.std(tin_array[processed_file,index,1]),\
                                        np.std(tin_array[processed_file,index,2]))])
        sample_median_tins.append([np.median(tin_array[processed_file,index,i]) for i in range(3)])                                
        sample_tin_sds.append([np.std(tin_array[processed_file,index,i]) for i in range(3)])
        sample_names.append(file.split('/')[-1][:-4])
        processed_file += 1
    
    SUM_SAMPLE.close()
    sample_median_tins = np.array(sample_median_tins)
    sample_tin_sds = np.array(sample_tin_sds)
    
    #Plot the median TINs of samples.
    index = np.arange(bamfile_number)
    if bamfile_number > 50:
        plt.figure(figsize = [bamfile_number/5, 20])
    else:
        plt.figure(figsize = [10,20])
    plt.subplot(411)   
    plt.bar(index, sample_median_tins[:,0], yerr=sample_tin_sds[:,0], color='deepskyblue')
    plt.xticks([])    
    plt.title('Median base-level TIN of samples')
    plt.subplot(412)
    plt.bar(index, sample_median_tins[:,1], yerr=sample_tin_sds[:,1], color='deepskyblue')
    plt.xticks([])
    plt.ylabel('nonuniform coverage <-- median TIN --> uniform coverage') 
    plt.title('Median exon-level TIN of samples')
    plt.subplot(413)
    plt.bar(index, sample_median_tins[:,2], yerr=sample_tin_sds[:,2], color='deepskyblue')
    plt.xticks([])   
    plt.title('Median exon-average TIN of samples')
    plt.subplot(414)
    plt.bar(index, expressed_transcript, color='deepskyblue')   
    plt.xticks(index, sample_names, rotation=90, fontsize=5)
    plt.ylabel('number')   
    plt.title('Number of expressed transcripts of samples')
    plt.savefig('summary_sample.pdf')
    plt.clf()
    
def summary_transcript():
    '''
    Summary KS and TIN values of each transcript in different samples.
    Output a transcript summary file and a plot pdf file in the current directory.
    '''
    #Summary file of high expression transcripts.
    SUM_TRANSCRIPT = open('summary_transcript.xls','w')
    print >>SUM_TRANSCRIPT, '\t'.join(['transcript','chrom',\
                                       'base_level_KS(mean)','exon_level_KS(mean)','exon_average_KS(mean)',\
                                       'base_level_TIN(median)','exon_level_TIN(median)','exon_average_TIN(median)',\
                                       'base_level_KS(sd)','exon_level_KS(sd)','exon_average_KS(sd)',\
                                       'base_level_TIN(sd)','exon_level_TIN(sd)','exon_average_TIN(sd)'])   
    processed_transcript = 0
    tins = []
    kss = [] 
    for line in open(options.ref_bed,'r'):
        if line.startswith(('#','track','browser')):continue  
        #Parse fields from BED file. 
        fields = line.split()
        chrom = fields[0]
        name = fields[3]
        #Only output transcript with high expression in more than certain ratio of samples.
        if np.mean(expression_level[:,processed_transcript]==2) >= options.minimum_sample_ratio:
            print >>SUM_TRANSCRIPT, '\t'.join( [str(i) for i in (name, chrom,\
                                                np.mean(ks_array[:,processed_transcript,0]),\
                                                np.mean(ks_array[:,processed_transcript,1]),\
                                                np.mean(ks_array[:,processed_transcript,2]),\
                                                np.median(tin_array[:,processed_transcript,0]),\
                                                np.median(tin_array[:,processed_transcript,1]),\
                                                np.median(tin_array[:,processed_transcript,2]),\
                                                np.std(ks_array[:,processed_transcript,0]),\
                                                np.std(ks_array[:,processed_transcript,1]),\
                                                np.std(ks_array[:,processed_transcript,2]),\
                                                np.std(tin_array[:,processed_transcript,0]),\
                                                np.std(tin_array[:,processed_transcript,1]),\
                                                np.std(tin_array[:,processed_transcript,2]))])
            tins.append([np.median(tin_array[:,processed_transcript,i]) for i in range(3)])
            kss.append([np.mean(ks_array[:,processed_transcript,i]) for i in range(3)])
        processed_transcript += 1
        
    SUM_TRANSCRIPT.close()
    kss = np.array(kss)
    tins = np.array(tins)
    
    #Plot the mean KS and median TIN of transcripts.
    plt.figure(figsize=[20,16])
    plt.subplot(221)
    plt.scatter(kss[:,0], tins[:,0], s=3)
    plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
    plt.xlabel("5' bias <-- KS --> 3' bias")
    plt.title('Mean base-level KS and median base-level TIN across samples')
    plt.subplot(223)
    plt.scatter(kss[:,1], tins[:,1], s=3)
    plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
    plt.xlabel("5' bias <-- KS --> 3' bias")
    plt.title('Mean exon-level KS and median exon-level TIN across samples')
    plt.subplot(224)
    plt.scatter(kss[:,2], tins[:,2], s=3)
    plt.ylabel('nonuniform coverage <-- TIN --> uniform coverage')
    plt.xlabel("5' bias <-- KS --> 3' bias")
    plt.title('Mean exon-average KS and median exon-average TIN across samples')
    plt.savefig('summary_transcript.pdf')
    plt.clf()

def summary_exon():
    '''
    Summary KS and TIN values of each transcript in different samples.
    Output an exon summary file in the current directory.
    '''
    #Summary file of high expression transcripts.
    SUM_EXON = open('summary_exon.xls','w')
    print >> SUM_EXON, '\t'.join(['exon','chrom','intra_exon_KS(mean)','intra_exon_TIN(median)',\
                                  'intra_exon_KS(sd)','intra_exon_TIN(sd)'])
    processed_exon = 0
    for line in open(options.ref_bed,'r'):
        if line.startswith(('#','track','browser')):continue  
        #Parse fields from BED file. 
        fields = line.split()
        chrom = fields[0]
        name = fields[3]
        exon_number = int(fields[9])
        for i in range(exon_number):
            #Only output exon expressed in more than certain ratio of samples.
            if np.mean(exon_array[:,processed_exon,0]!=0) >= options.minimum_sample_ratio:
                print >> SUM_EXON, '\t'.join([name+'.'+str(i+1), chrom,\
                                              str(np.mean(exon_array[:,processed_exon,0])),\
                                              str(np.median(exon_array[:,processed_exon,1])),\
                                              str(np.std(exon_array[:,processed_exon,0])),\
                                              str(np.std(exon_array[:,processed_exon,1]))])
            processed_exon += 1        
    SUM_EXON.close()
  
def output_rank():
    '''
    Output a file reporting TIN rank in all samples of transcripts. 
    '''
    for i in range(3):
        rank_array = np.apply_along_axis(rankdata, 0, tin_array[:,:,i])
        if i == 0:
            rank_file = 'base_level_TIN_rank.xls'
        elif i == 1:
            rank_file = 'exon_level_TIN_rank.xls'
        else:
            rank_file = 'exon_average_TIN_rank.xls'
        TIN_RANK = open(rank_file,'w')
        print >>TIN_RANK, '\t'.join(['gene'] + [j.split('/')[-1][:-4] for j in bamfiles])
        processed_transcript = 0
        for line in open(options.ref_bed,'r'):
            if line.startswith(('#','track','browser')):continue  
            #Parse fields from BED file. 
            fields = line.split()
            name = fields[3]
            if np.all(expression_level[:,processed_transcript] > 0):
                print >>TIN_RANK, '\t'.join([name] + [str(j) for j in rank_array[:,processed_transcript]])
            processed_transcript += 1
        TIN_RANK.close()
        

#Set options.
usage='%prog [options]' + '\n' + __doc__ + '\n'
parser = OptionParser(usage)
parser.add_option('-i','--input',action='store',type='string',dest='input_files',help="Input BAM file(s). Must be one of followings: 1) a single BAM file. 2) ',' separated BAM files (no spaces allowed). 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam files (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools. [required]")
parser.add_option('-r','--refbed',action='store',type='string',dest='ref_bed',help='Reference gene model in BED format. Must be strandard 12-column BED file. [required]')
parser.add_option('-d','--minDepth',action='store',type='float',dest='minimum_average_depth',default=0.5,help='Minimum average depth for each transcript to be processed. default=%default')
parser.add_option('-s','--minSample',action='store',type='float',dest='minimum_sample_ratio',default=1.0,help='Minimum ratio of samples expressing a transcript/exon. Only transcript/exon expressed in greater than this ratio of samples will be summarized. default=%default')
parser.add_option('-e','--exon',action='store_true',dest='exon',help='Output a xls file reporting intra-exon KS and TIN of each exon.')
parser.add_option('-k','--rank',action='store_true',dest='rank',help='Output a xls file reporting the rank of median transcript TIN across samples.')
(options,args)=parser.parse_args()

if not (options.input_files and options.ref_bed):
    parser.print_help()
    sys.exit(0)
    
if not os.path.exists(options.ref_bed):
    print >>sys.stderr, '\n\n' + options.ref_bed + ' does NOT exists' + '\n'
    parser.print_help()
    sys.exit(0)

#Read BAM files.    
printlog('Get BAM file ...')
bamfiles = getBamFiles.get_bam_files(options.input_files)
transcript_number, exon_number = count_number(options.ref_bed)
bamfile_number = len(bamfiles)
    
if bamfile_number <= 0:
    print >>sys.stderr, 'No BAM file found, exit.'
    sys.exit(0)
else:
    print >>sys.stderr, 'Total %d BAM file(s):' % len(bamfiles)
        
#Matrix of ks and tin of all transcripts in all samples.
ks_array = np.zeros([bamfile_number, transcript_number, 3])
tin_array = np.zeros([bamfile_number, transcript_number, 3])
#Matrix of expression level of all transcripts in all samples.
expression_level = np.zeros([bamfile_number, transcript_number])
if options.exon:
    #Matrix of ks and tin of all transcripts in all samples.
    exon_array = np.zeros([bamfile_number, exon_number, 2])
    
process_bamfiles(bamfiles)
summary_transcript()
if options.exon:
    summary_exon()
if options.rank:
    output_rank()