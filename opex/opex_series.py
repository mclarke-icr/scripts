#!/usr/bin/env python2.7

"""
Last updated: 30/09/2015
"""
import os
import argparse
from collections import OrderedDict
import pysam
import multiprocessing

class SingleJob(multiprocessing.Process):
    
    def __init__(self,sample,options,result_queue):
        multiprocessing.Process.__init__(self)
        
        #init variables
        self.sample = sample
        self.options = options
        self.result_queue = result_queue
        
        #read in files
        if len(self.options.bam) > 0:
            self.samfile=pysam.Samfile(self.sample+self.options.bam, "rb")  
        self.variants=[ line.rstrip() for line in open(self.sample+"_annotated_calls.txt","r") ]
        self.header =self.variants.pop(0)

    def getCoverage(self,chrom,pos):
        #allowing for both chr# and # noteation of chromosomes 
        
        #chrprefix=self.samfile.references[0].startswith('chr')
        #if chrprefix and not chrom.startswith('chr'): goodchrom='chr'+chrom
        #if not chrprefix and chrom.startswith('chr'): goodchrom=chrom[3:]
        #print "___"+pos+"___"
        puo = self.samfile.pileup(chrom,int(pos),int(pos)+1)
        print len(puo)
        for pileupcolumn in puo:
            #print pileupcolumn.pos
            if pileupcolumn.pos == pos:
                #print pos
                cov = pileupcolumn.n
                
                if cov >=15 :
                    return 1
                else:
                    return 0

    def run(self):
        var_summary = {}
        het_gt = ["0/1","1/0","1/2","2/1"]
        multi_gt = ["1/2","2/1"]
        for var in self.variants:
            summary_fields =[]
            fields = var.split("\t")
            #geting coverage above 15x
            cov15 = self.getCoverage(fields[0],fields[1])
            #variant specific index "gene_transcript_CSN"
            idx_id = "_".join([fields[13],fields[12],fields[16]])
            #first time variant is seen
            if idx_id not in var_summary:
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
                summary_fields.append(cov15)
                var_summary[idx_id] = summary_fields
            # subsequent times variant is seen
            else:
                summary_fields = var_summary[idx_id]
                if fields[10] == "1/1":
                    if fields[5] == "high":
                        summary_fields[16] += 1 #high hom
                    else:
                        summary_fields[17] += 1 #lo hom
                elif fields[10] in het_gt:
                    if fields[5] == "high":
                        summary_fields[18] += 1 #high het
                    else:
                        summary_fields[19] += 1 #lo het
                if fields[10] in multi_gt:
                    summary_fields[20] += 1 #multi allelic
                    
                summary_fields[21] += cov15
                var_summary[idx_id] = summary_fields


        self.result_queue.put(var_summary)

        
#############################################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and out puts summary info from your opex runs')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-b","--bam")
parser.add_argument("-q","--quality",action='store_true')
parser.add_argument("-c","--class",action='store_true')
parser.add_argument("-t","-threads",default=1)
options = parser.parse_args()
print options.input

# filter lists
q_list = ["low"]
c_list =["INT","3PU","5PU"]
#setting up jobs
samples = [ line.rstrip() for line in open(options.input,"r") ]
job_list = []
tables = multiprocessing.Queue()
#for sample in samples:
#
#    mp_obj = SingleJob(sample,options,tables)
#    
#    job_list.append(mp_obj)
#
mp_obj = SingleJob(samples[0],options,tables)
job_list.append(mp_obj)

#running jobs
for job in job_list: job.start()
#for job in job_list: job.join()

#collecting returned data
all_var = {}
for sample in samples:
    outdata = tables.get()
    #print outdata
        
        