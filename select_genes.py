#!/usr/bin/env python2.7

"""
Last updated: 30/09/2015
"""
import os
import argparse



##############################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and a list of genes and goves the varainst in thoes genes for thoes samples, with exac and ICR100 annoations')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-g","--genes",required=True)
parser.add_argument("-o","--output",required=True)
parser.add_argument("-c","--icr1000",action='store_true',default=False)
parser.add_argument("-e","--exac",action='store_true',default=False)
parser.add_argument("-s", "--suffix",default="_annotated_calls.txt",type=str)

options = parser.parse_args()
print options.input

samples ={}
for line in open(options.input,"r"):
     l = line.rstrip() 
     sampinfo = l.split("\t")
     samples[sampinfo[0]] = sampinfo[1:]
gene_list = [ line.rstrip() for line in open(options.genes,"r") ]
outlines = []
for sample in samples:
    with open(workingdir+"/"+sample+options.suffix,"r") as file:
        for line in file:
            l = line.rstrip()
            varcols = l.split("\t")
            varcols[9] = sample
            if varcols[13] in gene_list:
                varcols.extend(samples[sample])
                outlines.append(varcols)
                
#adding annotation if neended
if options.exac:
    exac= {}
    with open("/scratch/cancgene/gst/databases/20151113_ExAC.r0.3.nonTCGA.sites_cavaOutput_uniqENSTandCSN.txt","r") as exac_fh:
        header=exac_fh.readline()
        h  = header.split("\t")

        for line in exac_fh:
            l = line.rstrip()
            var = l.split("\t")
            var_id = "_".join([var[4],var[5],var[6]]) # gene transcript csn
            if int(var[14]) ==0 :
                exac[var_id] = 0.0
            else:
                exac[var_id] = str((float(var[12])/float(var[14]))*100.0)
                
if options.icr1000:
    icr1000={}
    with open("/scratch/cancgene/gst/databases/20151008_Control_unionOfAllVars_counts.txt","r") as icr1000_fh:
        header=icr1000_fh.readline()
        h  = header.split("\t")
        for line in icr1000_fh:
            l = line.rstrip()
            var = l.split("\t")
            var_id = "_".join([var[6],var[5],var[9]]) # gene transcript csn
            icr1000[var_id] = str((float(var[16])/993.0)*100.0)
            

out=open(options.output,"w")
for line in outlines:
    var_id = "_".join([line[13],line[12],line[16]])
    if options.exac:
        if var_id in exac:
            line.append(exac[var_id])
        else:
            line.append("0")
    if options.icr1000:
        if var_id in icr1000:
            line.append(icr1000[var_id])
        else:
            line.append("0")
    
    outline= "\t".join(line)+"\n"
    out.write(outline)

    
    