#!/usr/bin/env python2.7

"""
Last updated: 30/09/2015
"""
import os
import argparse
#from collections import OrderedDict
import pysam
from multiprocessing import Pool


def run(sample):

    var_summary = {}
    het_gt = ["0/1","1/0","1/2","2/1"]
    multi_gt = ["1/2","2/1"]
    variants=[ line.rstrip() for line in open(sample+"_annotated_calls.txt","r") ]
    header =variants.pop(0)
    for var in variants:
        summary_fields =[]
        fields = var.split("\t")
        #geting coverage above 15x
        if len(options.bam) > 0:
            cov15 = 0
            samfile=pysam.Samfile(sample+options.bam, "rb") 
            puo = samfile.pileup(fields[0],int(fields[1]),int(fields[1])+1)
            for pileupcolumn in puo:
                if pileupcolumn.pos == int(fields[1]):
                    cov = pileupcolumn.n
                    if cov >=15 : cov15 = 1
        #variant specific index "gene_transcript_CSN"
        idx_id = "_".join([fields[13],fields[12],fields[16]])
        summary_fields =[]
        # adding variant specific fields
        summary_fields.extend(fields[0:4])
        summary_fields.extend(fields[11:23])
        #adding counts
        countlist = [0,0,0,0,0]
        if fields[10] == "1/1":
            if fields[5] == "high":
                countlist[0] = 1 #high hom
            else:
                countlist[1] = 1 #lo hom
        elif fields[10] in het_gt:
            if fields[5] == "high":
                countlist[2] = 1 #high het
            else:
                countlist[3] = 1 #lo het
        if fields[10] in multi_gt:
            countlist[4] = 1 #multi allelic
        
        summary_fields.extend(countlist)
        
        if len(options.bam) > 0:
            summary_fields.append(cov15)
        else:
            summary_fields.append("N/A")
        #adding "TR,TC" to varinat specific list
        summary_fields.append(",".join([fields[7],fields[8]]))
        
        var_summary[idx_id] = summary_fields

    return var_summary
    
        
#############################################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and out puts summary info from your opex runs')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-b","--bam")
parser.add_argument("-q","--quality",action='store_true')
parser.add_argument("-c","--class",action='store_true')
parser.add_argument("-t","--threads",default=1,type=int)
options = parser.parse_args()
print options.input

# filter lists
q_list = ["low"]
c_list =["INT","3PU","5PU"]

#setting up jobs
samples = [ line.rstrip() for line in open(options.input,"r") ]

job = Pool(options.threads)

returned_list = job.map(run,samples)

#running jobs
"""
#collecting returned data
all_var = {}
samp_var_m = {}
nsamp = len(samples)
for samp_idx in range(len(samples)):
    outdata = tables.get()
    for idx in outdata:
        if idx not in all_var:
            all_var[idx] = outdata[idx] 
            # adding row to samp*var matrix
            samp_var_m[idx] = [0] * nsamp
            samp_var_m[idx][samp_idx] = outdata[idx][22]
        else:
            # summing number 
            all_var[idx][16] += outdata[idx][16]
            all_var[idx][17] += outdata[idx][17]
            all_var[idx][18] += outdata[idx][18]
            all_var[idx][19] += outdata[idx][19]
            all_var[idx][20] += outdata[idx][20]
            all_var[idx][21] += outdata[idx][21]
        
            samp_var_m[idx][samp_idx] = outdata[idx][22]
        
        
        
out = open( ".".join(options.input.split(".")[:-1])+"_variant_summary_individuals.txt", "w")
for idx in all_var:
    var_sum = "\t".join(all_var[idx][0:21])
    cov15 = round( all_var[idx][21] / nsamp *100 , 1)
    out_line = var_sum+"\t"+str(cov15)+"\n"
    
"""