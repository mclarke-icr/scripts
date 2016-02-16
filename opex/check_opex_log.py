#!/usr/bin/env python2.7

import sys

all_opex_log = sys.argv[1]
out=open("./failed.txt","w")
for line in open(all_opex_log,"r"):
    l = line.rstrip()
    cols = l.split("\t")
    times = [int(t)for t in cols[1:]]
    if len(times)== 8:
        if sum(times[1:]) < 10000:  
            out.write(str(cols[0])+" failed, finished too fast\n")
    elif len(times)== 7:
        out.write(str(cols[0])+" failed at Creating output txt file\n")
    elif len(times)== 6:
        out.write(str(cols[0])+" failed at Variant annotation\n")
    elif len(times)== 5:
        out.write(str(cols[0])+" failed at Variant calling\n")
    elif len(times)== 4:
        out.write(str(cols[0])+" failed at Checking coverage\n")
    elif len(times)== 3:
        out.write(str(cols[0])+" failed at Duplicate marking\n")
    elif len(times)== 2:
        out.write(str(cols[0])+" failed at Comparing number of reads in FASTQ and BAM files\n")
    elif len(times)== 1:
        out.write(str(cols[0])+" failed at Mapping reads, converting output to BAM file, sorting and indexing BAM\n")
    else:
        out.write(str(cols[0])+" failed\n")

