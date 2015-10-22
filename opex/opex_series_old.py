#!/usr/bin/env python2.7

"""
Last updated: 30/09/2015
"""
import os
import argparse
from collections import OrderedDict
import pysam
# delim object
class Delim(object):
    def __init__(self,path,hea,dlm):
        self.path = path
        file_matrix = [] 
        file_lines = []
        header = "NA"
        with open(path,"r") as f:
            header = f.readline()
            cols = header.split(dlm)
            if hea == "header":
                for col in range(len(cols)):
                    file_matrix.append([])
            else:
                hln = header.rstrip('\r\n')
                file_lines.append(hln)
                for c in range(len(cols)):
                    file_matrix.append([])
                    file_matrix[c].append(cols[c])
            for line in f:
                ln = line.rstrip('\r\n')
                file_lines.append(ln)
                fields = ln.split(dlm)  
                for i in range(len(fields)):
                    file_matrix[i].append(fields[i])
        self.data = file_matrix # table of the file .data[x][y] where x is the col and y is the row
        self.lines = file_lines
        self.header = header
        
    def get_gene_trans_csn(self):
        _gene_trans_csn=[]
        for i in range(len(self.data[0])):
            _gene_trans_csn.append("_".join([self.data[13][i],self.data[16][i]]))
        return _gene_trans_csn

#############################################################################################################
scriptdir=os.path.dirname(os.path.realpath(__file__))
workingdir=os.getcwd()

#command line args
parser = argparse.ArgumentParser(description='takes in a list of samples and out puts summary info from your opex runs')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-q","--quality",action='store_true')
parser.add_argument("-c","--class",action='store_true')
args = vars(parser.parse_args())
infile = args["input"]
qual_filt = args["quality"]
class_filt = args["class"]
# filter lists
q_list = ["low"]
c_list =["INT","3PU","5PU"]

print infile
if infile.startswith("/"):
    base_dir ="/"+"/".join(infile.split("/")[:-1])
else:
    base_dir = "/".join(infile.split("/")[:-1])
    base_dir=workingdir+"/"+base_dir

samples = [ line.rstrip() for line in open(infile,"r") ]

# initiating collumns
variants ={}
samp_m = OrderedDict()
n_samp = len(samples)

het_gt = ["0/1","1/0","1/2","2/1"]
multi_gt = ["1/2","2/1"]
for samp_idx in range(len(samples)):
    samp = Delim(base_dir+"/"+samples[samp_idx]+"_annotated_calls.txt","header","\t")

    for idx in range(len(samp.lines)):
        if class_filt == True and samp.data[17][idx] in c_list and samp.data[21][idx] in c_list:
            continue
        if qual_filt == True and samp.data[4][idx] in q_list:
            continue
        idx_id = "_".join([samp.data[13][idx],samp.data[12][idx],samp.data[16][idx]])

        if idx_id not in variants :
            # variant sample matrix
            #samp_m[idx_id] = [0] * n_samp
            #samp_m[idx_id][samp_idx] = "/".join([samp.data[7][idx],samp.data[8][idx]])
            
            varlist = [] # initiate list
            
            varlist.append(samp.data[0][idx]) #chr
            varlist.append(samp.data[1][idx]) #pos
            varlist.append(samp.data[2][idx]) #ref
            varlist.append(samp.data[3][idx]) #alt

            varlist.append(samp.data[11][idx]) #type
            varlist.append(samp.data[12][idx]) #transcript
            varlist.append(samp.data[13][idx]) #gene
            varlist.append(samp.data[14][idx]) #trinfo
            varlist.append(samp.data[15][idx]) #loc
            varlist.append(samp.data[16][idx]) #csn
            varlist.append(samp.data[17][idx]) #class
            varlist.append(samp.data[18][idx]) #so
            varlist.append(samp.data[19][idx]) #impact
            varlist.append(samp.data[20][idx]) #altann
            varlist.append(samp.data[21][idx]) #altclass
            varlist.append(samp.data[22][idx]) #altso
            
            countlist = [0,0,0,0,0]
            if samp.data[10][idx] == "1/1":
                if samp.data[5][idx] == "high":
                    countlist[0] = 1 #high hom
                else:
                    countlist[1] = 1 #lo hom
            elif samp.data[10][idx] in het_gt:
                if samp.data[5][idx] == "high":
                    countlist[2] = 1 #high het
                else:
                    countlist[3] = 1 #lo het
            if samp.data[10][idx] in multi_gt:
                countlist[4] = 1 #multi allelic
            
            variants[idx_id] = varlist+countlist

        else:
            # variant sample matrix
            #samp_m[idx_id][samp_idx] ="/".join([samp.data[7][idx],samp.data[8][idx]])

            varlist = variants[idx_id]
            
            if samp.data[10][idx] == "1/1":
                if samp.data[5][idx] == "high":
                    varlist[16] = varlist[16] +1
                else:
                    varlist[17] = varlist[17] +1
            elif samp.data[10][idx] in het_gt :
                if samp.data[5][idx] == "high":
                    varlist[18] = varlist[18] +1
                else:
                    varlist[19] = varlist[19] +1
            if samp.data[10][idx] in multi_gt:
                varlist[20] = varlist[20] +1
            
            variants[idx_id] = varlist
            
        # sample matrix
outfile = ".".join(infile.split(".")[:-1])
out = open(outfile+"_variant_summary_individuals.txt","w")
out.write("CHR	POS	REF	ALT	TYPE	ENST	GENE	TRINFO	LOC	CSN	CLASS	SO	IMPACT	ALTCSN	ALTCLASS	ALTSO	HIHOM	LOHOM	HIHET	LOHET	MULTIALLELE\n")

for idx_id in variants:
    out_list = [ str(v) for v in variants[idx_id] ]
    #for samp_idx in samp_m[idx_id]:
        #print samp_m[idx_id]
    out_line = "\t".join(out_list)+"\n"
    #print out_line
    out.write(out_line)
    
"""    
CHROM	0
POS	1
REF	2
ALT	3
QUAL	4
QUALFLAG	5
FILTER	6
TR	7
TC	8
SAMPLE	9
GT	10
TYPE	11
ENST	12
GENE	13
TRINFO	14
LOC	15
CSN	16
CLASS	17
SO	18
IMPACT	19
ALTANN	20
ALTCLASS	21
ALTSO	22

Opex_series.py
Input - list of samples for analysis. Assumes that each sample has been analysed with OpEx thus has expected suffixes, i.e. _annotated_calls.txt and _picard.bam/_picard.bai
Output
1)	_variant_summary.txt
    a.	Tab-separated
    b.	Rows
        i.	All variants seen in the series
        ii.	No duplicates as defined by Gene + ENST + CSN
    c.	Columns
        i.	OpEx generated
            1.	TYPE
            2.	ENST
            3.	GENE
            4.	TRINFO
            5.	LOC
            6.	CSN
            7.	CLASS
            8.	SO
            9.	IMPACT
            10.	ALTANN
            11.	ALTCLASS
            12.	ALTSO
        ii.	OpEx generated, selected from first occurrence
            1.	CHROM
            2.	POS
            3.	REF
            4.	ALT
        iii.	Counts
            1.	# high quality homozygous (GT=1/1, QUALFLAG="high") 
            2.	# high quality homozygous (GT=0/1 or 1/0, QUALFLAG="high") 
            3.	# low quality homozygous (GT=1/1, QUALFLAG="low") 
            4.	# low quality homozygous (GT=0/1 or 1/0, QUALFLAG="low") 
            5.	# multi allelic (GT is not 1/1, 0/1, or 1/0)
    d.	Optional filters
        i.	Exclude those with CLASS=INT/3PU/5PU and ALTCLASS=INT/3PU/5PU/"."
        ii.	Exclude variants where all calls are low quality
2)	_variant_summary_individuals.txt
    a.	Rows
        i.	Variants, exactly corresponding with those in _variant_summary.txt
    b.	Columns
        i.	The samples supplied by the user
        ii.	Each field is either 0 if the variant was not seen in the sample, or "TR/TC" from the sample's file for that variant.
3)	_variant_summary_denominators.txt - ER: this could be an option?? Maybe easier as entirely separate given BAM requirement
    a.	Rows
        i.	Variants, exactly corresponding with those in _variant_summary.txt
    b.	Column
        i.	# of samples with >=15x coverage, retrieved in the same way as denovo.py

"""       
        
        
        
        
        
        
        
        
        