#!/usr/bin/python
import re
import sys
sys.path.append("/scratch/cancgene/mclarke/python/lib")
import biofile
#command line args, in and out file names
infile = sys.argv[1]
family = sys.argv[2]
print infile
base_dir = "/".join(infile.split("/")[:-1])
print  base_dir
outfile = family+"_Denovo.txt"
print outfile
out = open(outfile,"w")
# extracting trio info
#"""
info = biofile.Delim(infile,"header","\t")

mother_re = re.compile("[Mm]other")
father_re = re.compile("[Ff]ather")
proband_re = re.compile("[Pp]roband")

out = open(outfile,"w")

out.write("\t".join(["CHR","POS","REF","ALT","FLAG","GENO","SAMPLE","GENE","TRANS","EX","HGVS","CLASS","V1","DBSNP","V2","PTV","TYPE"])+"\n")

proband_list = []
for family in set(COGinfo.data[6]):
    p_idx = [ i for i in range(len(COGinfo.lines)) if proband_re.match(info.data[7][i]) and COGinfo.data[6][i] == family]
    m_idx = [ i for i in range(len(COGinfo.lines)) if mother_re.match(info.data[7][i]) and COGinfo.data[6][i] == family]
    f_idx = [ i for i in range(len(COGinfo.lines)) if father_re.match(COGinfo.data[7][i]) and COGinfo.data[6][i] == family]
    if len(p_idx)==1 and len(m_idx)==1 and len(f_idx)==1:
        mother = biofile.Dbtxt(base_dir+info.data[0][m_idx[0]]+".txt","nohea")
        father = biofile.Dbtxt(base_dir+info.data[0][f_idx[0]]+".txt","nohea")
        proband = biofile.Dbtxt(base_dir+info.data[0][p_idx[0]]+".txt","nohea")
        p_gene_hgvs,var_idx = proband.get_gene_hgvs()
        m_gene_hgvs,m_idx = mother.get_gene_hgvs()
        f_gene_hgvs,f_idx = father.get_gene_hgvs()

        mf_gene_hgvs = set(m_gene_hgvs + f_gene_hgvs)
        dn_idx =[]
        for p in range(len(var_idx)):
            total_var += 1
            if p_gene_hgvs[p] not in mf_gene_hgvs:
                   dn_idx.append(var_idx[p])
        for dn in set(dn_idx):
            out.write(proband.lines[dn]+"\n")
        roband_list.append(COGinfo.data[0][p_idx[0]])
print total_var

pbl = open("./probandlist.txt","w")

for pl in proband_list:
    pbl.write(pl+"\n")
            