#!/usr/bin/python

import sys
import re
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

################################################################################
infile =  sys.argv[1]
exclude_re = re.compile("[Dd]uplicate|[Ee]xcluded")
samples = [ line.rstrip() for line in open(infile,"r") ]

info = Delim("/scratch/cancgene/mclarke/REF/Exome_Processing_Key_20150707.txt", "header","\t")


for sample in samples:

    for idx in range(len(info.lines)):
        if exclude_re.match(info.data[16][idx]):
            continue
        else:
            if sample == info.data[2][idx]:
                #print "/data/sutlt/DGE/GENSUSC/gst/"+info.data[21][idx]+"/BAM/"+info.data[14][idx]+"/"+sample+"*.ba*"
                print "/data/sutlt/DGE/GENSUSC/gst/"+info.data[21][idx]+"/FASTQ/"+info.data[14][idx]+"/"+sample+"*.fastq.gz"
                #print sample+"\t"+info.lines[idx]
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            