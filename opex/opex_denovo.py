#!/usr/bin/python

"""
Last updated: 26/3/2015
"""
import re
import argparse

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
        

#########################################################################################################################################

#command line args
parser = argparse.ArgumentParser(description='selects variants only seen in Probands from Trio data sets')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-f","--family")
parser.add_argument("-o","--output")
args = vars(parser.parse_args())
infile = args["input"]
family = args["family"]
outfile = args["output"]

print infile
base_dir = "/".join(infile.split("/")[:-1])+"/"
if outfile == None:
    outfile = "".join(infile.split(".txt"))+"_denovo.txt"
print outfile
out = open(outfile,"w")

# extracting trio info

info = Delim(infile,"noheader","\t")

mother_re = re.compile("[Mm]other")
father_re = re.compile("[Ff]ather")
proband_re = re.compile("[Pp]roband")

out = open(outfile,"w")

out.write("\t".join(["CHR","POS","REF","ALT","FLAG","GENO","SAMPLE","GENE","TRANS","EX","HGVS","CLASS","V1","DBSNP","V2","PTV","TYPE"])+"\n")

total_var =0
if family == None:
    family = set(info.data[1])
else:
    family = family.split(",")
for fam in family:
    p_idx = [ i for i in range(len(info.lines)) if proband_re.match(info.data[2][i]) and info.data[1][i] == fam]
    m_idx = [ i for i in range(len(info.lines)) if mother_re.match(info.data[2][i]) and info.data[1][i] == fam]
    f_idx = [ i for i in range(len(info.lines)) if father_re.match(info.data[2][i]) and info.data[1][i] == fam]
    if len(p_idx)==1 and len(m_idx)==1 and len(f_idx)==1:
        mother = Delim(base_dir+info.data[0][m_idx[0]]+"_annotated_calls.txt","header","\t")
        father = Delim(base_dir+info.data[0][f_idx[0]]+"_annotated_calls.txt","header","\t")
        proband = Delim(base_dir+info.data[0][p_idx[0]]+"_annotated_calls.txt","header","\t")
        p_gene_hgvs = proband.get_gene_hgvs()
        m_gene_hgvs = mother.get_gene_hgvs()
        f_gene_hgvs = father.get_gene_hgvs()

        mf_gene_hgvs = set(m_gene_hgvs + f_gene_hgvs)
        dn_idx =[]
        for p in range(len(p_gene_hgvs)):
            total_var += 1
            if p_gene_hgvs[p] not in mf_gene_hgvs:
                   dn_idx.append(p)
        for dn in set(dn_idx):
            out.write(proband.lines[dn]+"\n")

print total_var
