#!/usr/bin/python

"""
Last updated: 25/09/2015
"""

import re
import argparse
from subprocess import call

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
        
    def get_gene_hgvs(self):
        _gene_hgvs=[]
        for i in range(len(self.data[0])):
            _gene_hgvs.append("_".join([self.data[13][i],self.data[16][i]]))
        return _gene_hgvs
        
# family info
class FamilyInfo(object):
    def __init__(self,path):
        info = open(path,"r")
        info_header = info.readline()
        family_info={}
        for line in info.readlines():
            (sample,family,relationship) = line.split("\t")
            rel = relationship.rstrip()
            if family in family_info:
                family_info[family][rel] = sample
            else:
                family_info[family] = { rel:sample }
        info.close()
        self.data = family_info
 
#########################################################################################################################################

#command line args
parser = argparse.ArgumentParser(description='selects variants only seen in Probands from Trio data sets')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-f","--family")
args = vars(parser.parse_args())
infile = args["input"]
family = args["family"]


print infile
base_dir = "/".join(infile.split("/")[:-1])
if len(base_dir) == 0:
    base_dir="."

# extracting trio info
info = FamilyInfo(infile)

mother_re = re.compile("[Mm]other")
father_re = re.compile("[Ff]ather")
proband_re = re.compile("[Pp]roband")

total_var =0
if family == None:
    family = set(info.data)
else:
    family = family.split(",")

for fam in family:
    # initiating out files
    bed = open(base_dir+"/"+fam+".bed","w")
    out = open(base_dir+"/"+fam+"_denovo.txt","w")
    firstline="#CHROM\tPOS\tREF\tALT\tQUAL\tQUALFLAG\tFILTER\tTR\tTC\tSAMPLE\tGT\tTYPE\tENST\tGENE\tTRINFO\tLOC\tCSN\tCLASS\tSO\tIMPACT\tALTANN\tALTCLASS\tALTSO\tCOV\tMOTHERCOV\tFATHERCOV\n"
    out.write(firstline)
    # get samp IDs
    pro = info.data[fam]['Proband']
    mo = info.data[fam]['Mother']
    fa = info.data[fam]['Father']
    # read in .txt files
    mother = Delim(base_dir+"/"+mo+"_annotated_calls.txt","header","\t")
    father = Delim(base_dir+"/"+fa+"_annotated_calls.txt","header","\t")
    proband = Delim(base_dir+"/"+pro+"_annotated_calls.txt","header","\t")
    # get gene_hgvs
    p_gene_hgvs = proband.get_gene_hgvs()
    m_gene_hgvs = mother.get_gene_hgvs()
    f_gene_hgvs = father.get_gene_hgvs()
    mf_gene_hgvs = set(m_gene_hgvs + f_gene_hgvs)
    # selecting denovo
    dn_idx =[]
    for p in range(len(p_gene_hgvs)):
        total_var += 1
        if p_gene_hgvs[p] not in mf_gene_hgvs:
               dn_idx.append(p)
    # writing .txt and bed files
    dn_set = set(dn_idx)
    for dn in dn_set:
        #out.write(proband.lines[dn]+"\n")
        bed.write("\t".join([proband.data[0][dn],str(int(proband.data[1][dn])-1),proband.data[1][dn],fam])+"\n")
    bed.close()
    #out.close()
    #running cava via bash
    exit_id = call(["python2.7","../opex-v1.0.0/tools/CoverView-v1.1.0/CoverView.py","-i",base_dir+"/"+pro+"_picard.bam", "-b", base_dir+"/"+fam+".bed", "-o", pro+"_dn_coverage"] )
    if exit_id != 0:
        exit("CoverView failiure: "+exit_id)
    exit_id = call(["python2.7","../opex-v1.0.0/tools/CoverView-v1.1.0/CoverView.py","-i",base_dir+"/"+mo+"_picard.bam", "-b", base_dir+"/"+fam+".bed", "-o", mo+"_dn_coverage"] )
    if exit_id != 0:
        exit("CoverView failiure: "+exit_id)
    exit_id = call(["python2.7","../opex-v1.0.0/tools/CoverView-v1.1.0/CoverView.py","-i",base_dir+"/"+fa+"_picard.bam", "-b", base_dir+"/"+fam+".bed", "-o", fa+"_dn_coverage"] )
    if exit_id != 0:
        exit("CoverView failiure: "+exit_id)
    #reading in coverage info  
    cv_pro = Delim(pro+"_coverview_regions.txt","header","\t")
    cv_mo = Delim(mo+"_coverview_regions.txt","header","\t")
    cv_fa = Delim(fa+"_coverview_regions.txt","header","\t")
    #collecting coverage data
    cv_data = {}
    for idx in range(len(cv_pro.lines)):
        all_cov = "\t".join([cv_pro.data[6][idx],cv_mo.data[6][idx],cv_fa.data[6][idx]])
        if pro in cv_data:
            if cv_pro.data[1][idx] in cv_data[pro]:
                cv_data[pro][cv_pro.data[1][idx]][cv_pro.data[3][idx]] = all_cov  
            else:
                cv_data[pro][cv_pro.data[1][idx]] = { cv_pro.data[3][idx]:all_cov } 
        else:
            cv_data[pro] = { cv_pro.data[1][idx]:{ cv_pro.data[3][idx]:all_cov } }
    
    for dn in dn_set:
        out_line = "\t".join([proband.lines[dn].rstrip(), cv_data[pro]["chr"+proband.data[0][dn]][proband.data[1][dn]]])+"\n"
        out.write(out_line)
    out.close()
