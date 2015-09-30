#!/bin/env python

import argparse

### def ### 
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
        
class FamilyInfo(object):
    def __init__(self,path):
        info = open(path,"r")
        info_header = info.readline()
        _family_info={}
        for line in info.readlines():
            (sample,family,relationship) = line.split("\t")
            rel = relationship.rstrip()
            if family in _family_info:
                _family_info[family][rel] = sample
            else:
                _family_info[family] = { rel:sample }
        info.close()
        self.data = _family_info
        
### main ###
# inputs
parser = argparse.ArgumentParser(description='add_cov.py -i infile.txt -f familyinfo.txt -o outfile')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-o","--output")
parser.add_argument("-f","--family")
args = vars(parser.parse_args())
infile = args["input"]
outfile = args["output"]
infofile =args["family"]

# family info
family_info = FamilyInfo(infofile)
#variants

cv_data = {}
for family in family_info:
    pro = family_info[family]['Proband']
    mo = family_info[family]['Mother']
    fa = family_info[family]['Father']
    
    cv_pro = Delim(pro+"_coverview_regions.txt","header","\t")
    cv_mo = Delim(mo+"_coverview_regions.txt","header","\t")
    cv_fa = Delim(fa+"_coverview_regions.txt","header","\t")
    for idx in range(len(cv_pro.lines)):
        all_cov = "\t".join([cv_pro.data[6][idx],cv_mo.data[6][idx],cv_fa.data[6][idx]])
        if pro in cv_data:
            if cv_pro.data[1][idx] in cv_data[pro]:
                cv_data[pro][cv_pro.data[1][idx]][cv_pro.data[3][idx]] = all_cov  
            else:
                cv_data[pro][cv_pro.data[1][idx]] = { cv_pro.data[3][idx]:all_cov } 
        else:
            cv_data[pro] = { cv_pro.data[1][idx]:{ cv_pro.data[3][idx]:all_cov } }
            
if outfile == None:
    out = open("".join(infile.split(".txt"))+"_cov.txt","w")
else:
    out = open("".join(outfile.split(".txt"))+".txt","w")
with open(infile,"r") as f:
    infile_header = f.readline()
    out_header =infile_header.rstrip()+"\tCOV\tMOTHERCOV\tFATHERCOV\n"
    out.write(out_header)
    for line in f.readlines():
        var = line.split("\t")
        l =line.rstrip()
        out.write( l+"\t"+cv_data[var[9]]["chr"+var[0]][var[1]]+"\n")







