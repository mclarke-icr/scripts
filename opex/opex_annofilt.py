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

#tmplog = open("tmp.log","w")


#############################################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and out puts summary info from your opex runs')
parser.add_argument("-i","--input",required=True)
#parser.add_argument("-e","--exac",default=False)
#parser.add_argument("-c","--icr1000",default=False)
parser.add_argument("-s","--suffix",default='_annotated_calls.txt',type=str)
options = parser.parse_args()
print options.input

# intite some things
vclass3 = ["SS","SY","INT","5PU","3PU","."]
fam_samp_ids = [ line.rstrip() for line in open(options.input,"r") ]


options.icr1000 = "/scratch/cancgene/gst/databases/20151008_Control_unionOfAllVars_counts.txt"
options.exac = "/scratch/cancgene/gst/databases/20151113_ExAC.r0.3.nonTCGA.sites_cavaOutput_uniqENSTandCSN.txt"
options.sanger = "/scratch/cancgene/mclarke/REF/ProbandSangerAnnotation_20151203_MC.txt"
#exac
exac= {}
with open(options.exac,"r") as exac_fh:
    header=exac_fh.readline()
    h  = header.split("\t")

    for line in exac_fh:
        l = line.rstrip()
        var = l.split("\t")
        var_id = "_".join([var[4],var[5],var[6]]) # gene transcript csn
        if int(var[14]) ==0 :
            exac[var_id] = 0.0
        else:
            exac[var_id] = (float(var[12])/float(var[14]))*100.0
#icr1000
icr1000={}
with open(options.icr1000,"r") as icr1000_fh:
    header=icr1000_fh.readline()
    h  = header.split("\t")
    for line in icr1000_fh:
        l = line.rstrip()
        var = l.split("\t")
        var_id = "_".join([var[6],var[5],var[9]]) # gene transcript csn
        icr1000[var_id] = (float(var[16])/993.0)*100.0
#sanger
sanger={}
for line in open(options.sanger,"r"):
    l = line.rstrip()
    var = l.split("\t")
    if fam_samp_ids[0].startswith("COG"):
        var_id = "_".join([var[0],var[2],var[3]]) # COGid gene  csn
    else:
        var_id = "_".join([var[1],var[2],var[3]]) # sampleID gene  csn
    if var_id in sanger:
        sanger[var_id] = sanger[var_id]+","+var[4]
    else:
        sanger[var_id] = var[4]


vclass3 = ["SS","SY","INT","5PU","3PU","."]
fam_samp_ids = [ line.rstrip() for line in open(options.input,"r") ]
head_todo = True
output=".".join(options.input.split(".")[:-1])+"_annofilt.txt"
out =open(output,"w")
for id in fam_samp_ids:
    fam_samp,pheno = id.split("\t")
    for line in open(fam_samp+options.suffix,"r") :
        outline =line.rstrip()
        var = outline.split("\t")
        #writing header
        if outline.startswith('#'):
            if head_todo:
                out.write(outline+"\tEXAC\tICR1000\tSANGER\tPHENO\n")
                head_todo =False
            continue
        var[9] = fam_samp
        outline = "\t".join(var)
        # filtering out class 3 variants
        if var[17] in vclass3 and var[21] in vclass3:
            continue
        var_id = "_".join([var[13],var[12],var[16]]) # gene transcript csn
        sanger_id = "_".join([var[9],var[13],var[16]]) # gene csn
        #filtering on allele freq in exac
        if var_id in exac:
            if exac[var_id] <= 0.5:
                outline = outline+"\t"+str(exac[var_id])
            else:
                continue
        else:
            outline = outline+"\t0"
        #filtering on allele freq in icr1000
        if var_id in icr1000 :
            if icr1000[var_id] <= 0.5:
                outline = outline+"\t"+str(icr1000[var_id])
            else:
                continue
        else:
            outline = outline+"\t0"
        #sanger confirmation annotation
        if sanger_id in sanger:
            outline = outline+"\t"+str(sanger[sanger_id])
        else:
            outline = outline+"\t."
        
        out.write(outline+"\t"+pheno+"\n")
        
    