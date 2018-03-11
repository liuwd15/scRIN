#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
This program calculates KS, TIN and some other metrics for each transcript or exon in BED file.

Copyright (C) 2018 Wendao Liu
Copyright (C) 2016 original author Liguo Wang

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
-------------------------------------------------------------------------------------------------'''
import sys,os
import math
import numpy
import pysam
from optparse import OptionParser
from qcmodule import getBamFiles
from time import strftime
from scipy.special import smirnov


def printlog (message):
    '''
    message: string.
    Print message into stderr with time string appending to it.
    '''
    message="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + message
    print >>sys.stderr,message

def genomic_positions(refbed,exon):
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
            chromStart = int( fields[1] )
            chromEnd = int( fields[2] )
            name = fields[3]
            strand = fields[5]
            exonCount = int(fields[9])               
            exonStart = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
            exonStart = map((lambda x: x + chromStart ), exonStart)
            exonEnd = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
            exonEnd = map((lambda x, y: x + y ), exonStart, exonEnd)
        except:
            print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
            continue
            
        #Return genomic positions of each nucleotide in exons.
        if exon:
            for i in range(exonCount):
                yield (name + '.' + str(i+1), chrom, exonStart[i], exonEnd[i], list(range(exonStart[i], exonEnd[i])), strand)
        #Return genomic positions of each nucleotide in transcripts.
        else:
            nucleotide_positions=[]
            for start,end in zip(exonStart,exonEnd):
                nucleotide_positions.extend(range(start,end))        
            yield (name, chrom, chromStart, chromEnd, nucleotide_positions, strand)
        

def calculate_coverage(samfile, chrom, positions):
    '''
    samfile: pysam.Samfile, input BAM file.
    chrom: string, index of chromosome.
    positions: list of integers, genomic positions of nucleotides.
    Return coverage on positions, a list of integers.
           coverage rate, a float number.
           peak, an integer.
    '''
    start = positions[0]            #Starting position of transcript/exon.
    end = positions[-1]             #Ending position of transcript/exon.
    length = len(positions)         #Length of transcript/exon.
    coverage = [0] * length         #Coverage of transcript/exon.
    covered_nucleotide = 0          #Number of covered nucleotides.
    peak = 0                        #Number of read peaks.
    i = 0
    
    #Calculate coverage of transcript/exon.
    for pileupcolumn in samfile.pileup(chrom, start, end):
        try:
            i = positions.index(pileupcolumn.pos,i)
            coverage[i] = pileupcolumn.n
            covered_nucleotide += 1
        except:
            continue
        
    coverage_rate = float(covered_nucleotide) / length       #Percentage of covered nucleotides.
    #Calculate number of read peaks.
    if coverage[0] > 0:
        peak += 1
    for i in range(1,length):
        if coverage[i] > 0 and coverage[i-1] == 0:
            peak += 1
    return coverage,coverage_rate,peak
    
def calculate_ks(coverage,read_length):
    '''
    coverage: list of integers, coverage on nucleotides.
    read_length: integer, the length of reads.
    Return KS and corresponding P value.
    '''
    length = len(coverage)
    depth = sum(coverage)            #depth = read number * read length
    #if length == 0 or depth == 0:
    #    return 0,1
    
    max_difference = 0.0 
    for i in range(1,length):
        #coverage here becomes empirical CDF of coverage.
        coverage[i] += coverage[i-1]
        if abs(max_difference) < abs(float(coverage[i])/depth - float(i+1)/length):
            max_difference = float(coverage[i])/depth - float(i+1)/length
            
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
    
    ks = math.sqrt(depth) * max_difference
    p = smirnov(int(depth/read_length),abs(max_difference))
    return ks,p

def calculate_tin(coverage):
    '''
    coverage: list of integers, coverage on nucleotides.
    Return TIN value.
    '''
    length = len(coverage)
    depth = sum(coverage)
    #if length == 0 or depth == 0:
    #    return 0
        
    entropy = 0.0
    for i in coverage:
        if i > 0:
            entropy += (float(i)/depth) * math.log(float(i)/depth)    
    tin = 100*(math.exp(-entropy)) / length
    return tin

#Set options.
usage="%prog [options]" + '\n' + __doc__ + '\n'
parser = OptionParser(usage)
parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input BAM file(s). "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files (no spaces allowed). 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam files (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools. [required]')
parser.add_option("-r","--refbed",action="store",type="string",dest="ref_bed",help='Reference gene model in BED format. Must be strandard 12-column BED file. [required]')
parser.add_option("-e","--exon",action="store_true",dest="exon",help="Calculate KS and TIN on exons rather than transcripts.")
parser.add_option("-d","--minDepth",action="store",type="int",dest="minimum_average_depth",default=1,help="Minimum average depth on each transcript/exon. default=%default")
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
    
if len(bamfiles) <= 0:
    print >>sys.stderr, "No BAM file found, exit."
    sys.exit(0)
else:
    print >>sys.stderr, "Total %d BAM file(s):" % len(bamfiles)
        
#Process each BAM file.
for f in bamfiles:
    print >>sys.stderr, "\t" + f         
    printlog("Processing " + f)
    samfile = pysam.Samfile(f, "rb")

    #Summary file of a sample.
    SUM = open(os.path.basename(f).replace('bam','') + 'summary.txt','w')
    print >>SUM, "\t".join(['Bam_file','ks(mean)','ks(stdev)','tin(median)','tin(stdev)'])
    
    #Details of each transcript/exon.
    OUT = open(os.path.basename(f).replace('bam','') + 'metrics.txt','w')
    print >>OUT, "\t".join(["geneID","chrom", "length","average_depth","coverage_rate","peak","tin","ks","p"])
    
    all_ks = []
    all_tin = []    
    finish = 0            #number of processed transcript/exon.
    read_length = samfile.head(1).next().infer_query_length()
    gene_object = "exons" if options.exon else "transcripts"
    
    #Process each transcript/exon.
    for name, chrom, chromStart, chromEnd, positions, strand in genomic_positions(options.ref_bed,options.exon):
        finish += 1
        length = len(positions) 
        coverage,rate,peak = calculate_coverage(samfile, chrom,sorted(positions))
        
        #Omit transcripts/exons with average depth.
        if numpy.mean(coverage) < options.minimum_average_depth:
            print >>OUT, '\t'.join([str(i) for i in (name, chrom, length, numpy.mean(coverage), rate, peak, 0, 0, 1)])
        #Process transcripts/exons above minimum average depth.
        else:
            ks,p = calculate_ks(coverage,read_length)
            tin = calculate_tin(coverage)
            #This step make sure that 3' bias results in positive KS and 5' bias results in negative KS.
            if strand == '+':
                ks = -ks
            all_ks.append(ks)
            all_tin.append(tin)
            print >>OUT, '\t'.join([str(i) for i in (name, chrom, length, numpy.mean(coverage), rate, peak, tin, ks,  p)])
        print >>sys.stderr, " %d %s finished\r" % (finish, gene_object),
    
    print >>SUM, "\t".join( [str(i) for i in (os.path.basename(f), numpy.mean(all_ks), numpy.std(all_ks), numpy.median(all_tin), numpy.std(all_tin))])
    OUT.close()
    SUM.close()
    samfile.close()
