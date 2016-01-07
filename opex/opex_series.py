#!/usr/bin/env python2.7

"""
Last updated: 30/09/2015
"""
import os
import argparse
#from collections import OrderedDict
import pysam
import multiprocessing
import time

tmplog = open("tmp.log","w")

class SingleJob(multiprocessing.Process):
    
    def __init__(self,fpath,options):
        multiprocessing.Process.__init__(self)
        
        #init variables
        self.fpath = fpath
        self.options = options
        #self.result_queue = result_queue

        #read in files
        #if self.options.bam != None:
        #    self.samfile=pysam.Samfile(self.fpath+self.options.bam, "rb")  
        self.variants=[ line.rstrip() for line in open(self.fpath,"r") ]
        self.header =self.variants.pop(0)

    def getCoverage(self,chrom,pos):
        #allowing for both chr# and # noteation of chromosomes 

        puo = self.samfile.pileup(chrom,int(pos),int(pos)+1)
        for pileupcolumn in puo:

            if pileupcolumn.pos == int(pos):
                #mapping and base quality for highquality
                #for pileupreadin in pileupcolumn.pileups:
                #    bq=pileupread.alignment.query_qualities[pileupread.query_position]
                #    mq=pileupread.alignment.mapq
                 
                cov = pileupcolumn.n
                
                if cov >=15 :
                    return 1
                else:
                    return 0

    def run(self):
        tmp_out = open(self.fpath+"_summary_tmp.txt","w")
        var_summary = {}
        het_gt = ["0/1","1/0","1/2","2/1"]
        multi_gt = ["1/2","2/1"]
        
        for var in self.variants:
            summary_fields =[]
            fields = var.split("\t")
            #geting coverage above 15x
            #if self.options.bam != None:
            #    cov15 = self.getCoverage(fields[0],fields[1])

            #variant specific index "gene_transcript_CSN"
            idx_id = "_".join([fields[13],fields[12],fields[16]])
            summary_fields =[]
            # adding variant specific fields
            summary_fields.append(idx_id)
            summary_fields.extend(fields[0:4])
            summary_fields.extend(fields[11:23])

            #adding counts
            countlist = ["0","0","0","0","0"]
            if fields[10] == "1/1":
                if fields[5] == "high":
                    countlist[0] = "1" #high hom
                else:
                    countlist[1] = "1" #lo hom
            elif fields[10] in het_gt:
                if fields[5] == "high":
                    countlist[2] = "1" #high het
                else:
                    countlist[3] = "1" #lo het
            if fields[10] in multi_gt:
                countlist[4] = "1" #multi allelic
            
            summary_fields.extend(countlist)

            #if self.options.bam != None and cov15 != None :
            #    summary_fields.append(str(cov15))
            #elif self.options.bam != None:
            #    summary_fields.append("0")
            #else:
            #    summary_fields.append("N/A")
            #adding "TR,TC" to varinat specific list
            summary_fields.append(",".join([fields[7],fields[8]]))

            outline = "\t".join(summary_fields)+"\n"
            tmp_out.write(outline)   

        return 0
#############################################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and out puts summary info from your opex runs')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-o","--output",required=True)
#parser.add_argument("-b","--bam",default=None)
#parser.add_argument("-q","--quality",action='store_true',default=False)
#parser.add_argument("-c","--class",action='store_true',default=False)
parser.add_argument("-t","--threads",default=1,type=int)
#parser.add_argument("-s", "--suffix",default="_annotated_calls.txt",type=str)
#parser.add_argument("-q","--basequality",default=10,type=int)
#parser.add_argument("-m","--mappingquality",default=20,type=int)
#parser.add_argument("-d","--incduplicates",action='store_true',default=False)
options = parser.parse_args()
print options.input

# filter lists
#q_list = ["low"]
#c_list =["INT","3PU","5PU"]

#setting up jobs
fpaths = [ line.rstrip() for line in open(options.input,"r") ]
job_list = []
#tables = multiprocessing.Queue()
print "make jobs"
for path in fpaths:

    mp_obj = SingleJob(path,options)
    
    job_list.append(mp_obj)

print "jobs made"
#running jobs
job_idx =0
running = 0
while job_idx < len(job_list):
    if running < options.threads:
        job_list[job_idx].start()
        job_idx += 1
    else:
        time.sleep(5)
    time.sleep(1)
    running = 0
    for job in job_list: 
        if job.is_alive(): running += 1
        
while len(multiprocessing.active_children()) > 0:
    time.sleep(10)
print "jobs run"
#collecting returned data
all_var = {}
samp_var_m = {}
nsamp = len(fpaths)
for samp_idx in range(len(fpaths)):
    tmpdata = open(fpaths[samp_idx]+"_summary_tmp.txt","r")
    for line in tmpdata:
        l = line.rstrip()
        outdata = l.split("\t")
        if outdata[22] == None:
            tmplog.write(line)

        if outdata[0] not in all_var:
            od_list = outdata[1:17]
            od_list.extend(map(int, outdata[17:22]))
            od_list.extend(outdata[22:23])
            all_var[outdata[0]] = od_list
            # adding row to samp*var matrix
            #samp_var_m[outdata[0]] = [0] * nsamp
            #samp_var_m[outdata[0]][samp_idx] = outdata[23]
        else:
            # summing number 
            all_var[outdata[0]][16] += int(outdata[17])
            all_var[outdata[0]][17] += int(outdata[18])
            all_var[outdata[0]][18] += int(outdata[19])
            all_var[outdata[0]][19] += int(outdata[20])
            all_var[outdata[0]][20] += int(outdata[21])
            #if options.bam != None :
            #    all_var[outdata[0]][21] += int(outdata[22])
            
            #samp_var_m[outdata[0]][samp_idx] = outdata[23]

        
out = open(options.output, "w")
for idx in all_var:
    vs =[ str(av) for av in all_var[idx][0:21] ]
    var_sum = "\t".join(vs)
    #if options.bam != None:
    #    cov15 = round( int(all_var[idx][21]) / nsamp *100 , 1)
    #else:
    #    cov15 ="N/A"
    #out_line = var_sum+"\t"+str(cov15)+"\n"
    out_line = var_sum+"\n"
    out.write(out_line)
out.close()
# individual matrix
#out = open( ".".join(options.input.split(".")[:-1])+"_variant_summary_individuals.txt", "w")
#out.write("var\t"+"\t".join(samples)+"\n")
#for idx in samp_var_m :
#    svm =[ str(samp_var_m[svl]) for svl in samp_var_m ]
#    out_line = idx+"\t"+"\t".join(str(samp_var_m[idx]))+"\n"
#    out.write(out_line)
#out.close()
for path in fpaths:
    os.remove(path+"_summary_tmp.txt")