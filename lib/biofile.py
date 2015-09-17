#!/usr/bin/python
"""
last modified:
2015-06-18

"""
#path ="/scratch/cancgene/mclarke/REF/Exome_Processing_Key_*.txt"
def mostRecent(path):
    import glob
    import re
    files = glob.glob(path)

    m = [ re.match("[0-9]{8}",f) for f in files ]
    return max([g.group() for g in m ],key=int)
    
class GlobalList(object):
    def __init__(self):
        self.list =[]
    def get_list_idx(self,new_list):
        list_set = set(self.list)
        out_idx = []
        for ls in new_list:
            if ls in list_set:
                out_idx.append(self.list.index( ls ))
            else:
                self.list.append( ls )
                out_idx.append(len(self.list)-1)
        return out_idx
"""     
# obsolete use GlobalList
class AllVar(object):
    def __init__(self):
        self.GH =[]
        self.CPRA =[]
        self.samp =[]
    def get_GH_idx(self,new_GH):
        GH_set = set(self.GH)
        out_idx = []
        for gh in new_GH:
            if gh in GH_set:
                out_idx.append(self.GH.index( gh ))
            else:
                self.GH.append( gh )
                out_idx.append(len(self.GH)-1)
        return out_idx

    def get_CPRA_idx(self,new_CPRA):
        CPRA_set = set(self.CPRA)
        out_idx = []
        for CPRA in new_CPRA:
            if CPRA in CPRA_set:
                out_idx.append(self.CPRA.index( CPRA ))
            else:
                self.CPRA.append( CPRA )
                out_idx.append(len(self.CPRA)-1)
        return out_idx
    
    def get_sample_idx(self,new_samp):
        samp_set = set(self.samp)
        out_idx = []
        for samp in new_samp:
            if samp in samp_set:
                out_idx.append(self.samp.index( samp ))
            else:
                self.samp.append( samp )
                out_idx.append(len(self.samp)-1)
        return out_idx
"""        
# The base  for delimited files
#requires: file path, the number of columns (or "header" if one is present), delimiter
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
     
    def get_cpra_custom(self,c_idx,p_idx,r_idx,a_idx,):
        _chr_pos_ref_alt=[]
        _idx=[]
        for i in range(len(self.data[c_idx])):
            r = self.data[r_idx][i].split(",")
            a = self.data[a_idx][i].split(",")
            if len(r) > len (a):
                for j in range(len(r)):
                    _chr_pos_ref_alt.append("_".join([self.data[c_idx][i],str(self.data[p_idx][i]),r[j],a[0]]))
                    _idx.append(i)
            elif len(r) < len (a):
                for j in range(len(a)):
                    _chr_pos_ref_alt.append("_".join([self.data[c_idx][i],str(self.data[p_idx][i]),r[0],a[j]]))
                    _idx.append(i)
            else:
                _chr_pos_ref_alt.append("_".join([self.data[c_idx][i],str(self.data[p_idx][i]),r[0],a[0]]))
                _idx.append(i)
        return _chr_pos_ref_alt,_idx

    def get_gh_custom(self,g_idx,h_idx):
        _gene_hgvs=[]
        _idx=[]
        for i in range(len(self.data[0])):
            g = self.data[g_idx][i].split(":")
            h = self.data[h_idx][i].split(":")
            for j in range(len(h)):
                ha = h[j].split(",")
                for k in range(len(ha)):
                    _gene_hgvs.append("_".join([g[j],ha[k]]))
                    _idx.append(i)
        return _gene_hgvs,_idx
        

class VCF(object):
    def __init__(self,path):
        import re
        self.head1_re = re.compile("##")
        self.head2_re = re.compile("#")
        self.path = path
        file_matrix = [] 
        file_lines = []
        header = ""
        with open(path,"r") as f:
            header = f.readline()
            while self.head1_re.match(header):
                header = f.readline()
            cols = header.split("\t")
            if self.head2_re.match(cols[0]):
                for col in range(len(cols)):
                    file_matrix.append([])
            else:
                for col in cols:
                    file_matrix.append([[col]])
            for line in f:
                file_lines.append(line)
                fields = line.split("\t")  
                for i in range(len(fields)):
                    file_matrix[i].append(fields[i])
        self.data = file_matrix # table of the file .data[x][y] where x is the col and y is the row
        self.lines = file_lines
        self.header = header

    def get_chr_pos_ref_alt(self):
        _chr_pos_ref_alt=[]
        _idx=[]
        for i in range(len(self.data[0])):
            r = self.data[3][i].split(",")
            a = self.data[4][i].split(",")
            if len(r) > len (a):
                for j in range(len(r)):
                    _chr_pos_ref_alt.append("_".join([self.data[0][i],self.data[1][i],r[j],a[0]]))
                    _idx.append(i)
            elif len(r) < len (a):
                for j in range(len(a)):
                    _chr_pos_ref_alt.append("_".join([self.data[0][i],self.data[1][i],r[0],a[j]]))
                    _idx.append(i)
            else:
                _chr_pos_ref_alt.append("_".join([self.data[0][i],self.data[1][i],r[0],a[0]]))
                _idx.append(i)
        return _chr_pos_ref_alt,_idx
        
    def list_info(self):
        _info_out =[]
        for info in self.data[7]:
            kvs = info.split(";")
            for kv in kvs:
                (k,v) = kv.split("=")
                _info_out.append(k)
        return set(_info_out)
    
    def get_info(self,key):
        import re
        _info_out =[]
        kvp_re = re.compile(";"+key+"=([^;]*)") 
        kvp_first_re = re.compile(key+"=([^;]*)")
        for info in self.data[7]:
            info_m = kvp_re.search(info)
            if info_m == None:
                info_m = kvp_first_re.match(info)
            _info_out.append(info_m.group(1))
        return _info_out
        


class Info(Delim):
    def __init__(self,path,hea):
        super(Info, self).__init__(path,hea,"\t")
        self.sample = self.data[0]
        self.mapped = self.data[1]
        self.cov15 = self.data[2]
        self.cov30 = self.data[3]
        self.cov50 = self.data[4]
        self.phenotype = self.data[5]
        self.family = self.data[6]
        self.relationship = self.data[7]
        self.fc = self.data[8]
        self.folder =self.data[9]
        self.trio = self.data[10]
        self.exomekit = self.data[11]
        
class Dbtxt(Delim):
    def __init__(self,path,hea):
        super(Dbtxt, self).__init__(path,hea,"\t")
        self.chr = self.data[0]
        self.pos = [int(x) for x in self.data[1]]
        self.ref = self.data[2]
        self.alt = self.data[3]
        self.quality = self.data[4]
        self.geno = self.data[5]
        self.sample = self.data[6]
        self.gene = self.data[7]
        self.ensembl_trans = self.data[8]
        self.exon = self.data[9]
        self.hgvs = self.data[10]
        self.var_class = self.data[11]
        self.alt_anno = self.data[12]
        self.db_snp = self.data[13]
        self.unk = self.data[14]
        self.ptv = self.data[15]
        self.var_type = self.data[16]
    # out puts orderd lists of each gene and hgvs from single and multi genic and allelic records also give a list of the original index.
    
    def get_gene_hgvs(self):
        _gene_hgvs=[]
        _idx=[]
        for i in range(len(self.chr)):
            g = self.gene[i].split(":")
            h = self.hgvs[i].split(":")
            for j in range(len(h)):
                ha = h[j].split(",")
                for k in range(len(ha)):
                    _gene_hgvs.append("_".join([g[j],ha[k]]))
                    _idx.append(i)
        return _gene_hgvs,_idx
        
    def get_chr_pos_ref_alt(self):
        _chr_pos_ref_alt=[]
        _idx=[]
        for i in range(len(self.chr)):
            r = self.ref[i].split(",")
            a = self.alt[i].split(",")
            if len(r) > len (a):
                for j in range(len(r)):
                    _chr_pos_ref_alt.append("_".join([self.chr[i],str(self.pos[i]),r[j],a[0]]))
                    _idx.append(i)
            elif len(r) < len (a):
                for j in range(len(a)):
                    _chr_pos_ref_alt.append("_".join([self.chr[i],str(self.pos[i]),r[0],a[j]]))
                    _idx.append(i)
            else:
                _chr_pos_ref_alt.append("_".join([self.chr[i],str(self.pos[i]),r[0],a[0]]))
                _idx.append(i)
        return _chr_pos_ref_alt,_idx
        
class Bed(Delim):
    def __init__(self,path,hea):
        super(Bed, self).__init__(path,hea," ")
        self.chr = self.data[0]
        self.start = [int(x) for x in self.data[1]]
        self.end = [int(x) for x in self.data[2]]
        
    #outputs list of position to which this applies, in vcf genomic coordinates
    #returns list of chr and pos lists
    def vcf_all_pos(self):
        chr_pos_list = []
        for n in range(0,len(self.chr)):
            for p in range(self.start[n]+1, self.end[n]):
                chr_pos_list.append([self.chr[n],p])

        return chr_pos_list
        
    #returns list of chr and pos lists in vcf genomic coordinates
    def pos_to_vcf(self):
        self.start = [ x+1 for x in self.start ]

        
class Tabix(object):
    def __init__(self,path):
        import sys
        sys.path.append("/users/cancgene/mclarke/.local/lib/python2.7/site-packages")
        import tabix
        self.tabix_obj = tabix.open(path)
    def get_range(self,chr,start,end):
        out_range=[]
        if start == end:
            tbx_range = self.tabix_obj.query(str(chr),start-1,end)
            is_first = True
            #print "___"
            for t in tbx_range:
                if start == int(t[1]):
                    out_range = [t]
        else:
            out_range = self.tabix_obj.query(str(chr),start,end)
        return out_range

class Pileup(Tabix):
    """
    reads in tabix indexed mpileup files
    the 'get_###' functions return specific columns between start and end positions
    any base without a record is assumed to have 0 coverage and is returned accordingly
    """
    def __init__(self,path):
        super(Pileup, self).__init__(path)
    """ ###broken###
    def get_all(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        prev_pos = start 
        for line in _tabix_out:
            while prev_pos < int(line[1]) -1:
                blank = [line[0],prev_pos,"-",0,0,0,0,0,""]
                out_list.append(blank)
                prev_pos += 1
            else:
                out_list.append(line)
                prev_pos += 1
        return out_list
    """
    def get_alt(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        cov = []
        pos = []
        for line in _tabix_out:
            cov.append(line[3])
            pos.append(line[1])
        idx =0
        for curr_pos in range(start,end+1):
            # if tabix outputs nothing
            if len(pos) ==0:
                out_list.append(["0","0","0","0"])
            # for covered positions
            elif int(pos[idx]) == curr_pos:
                out_list.append([line[4],line[5],line[6],line[7]])
                idx += 1
            # non-covered positions
            else:
                out_list.append(["0","0","0","0"])
                idx += 1
        return out_list
        
    def get_cov(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        cov = []
        pos = []
        for line in _tabix_out:
            cov.append(line[3])
            pos.append(line[1])
            
        idx =0
        for curr_pos in range(start,end+1):
            if len(pos) ==0:
                out_list.append('0')
            elif int(pos[idx]) == curr_pos:
                out_list.append(cov[idx])
                idx += 1
            else:
                out_list.append('0')
                idx += 1

        return out_list
    """ ###broken### fix using get_cov as template
    def get_a(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        prev_pos = start -1
        for line in _tabix_out:
            while prev_pos < int(line[1]) -1:
                out_list.append(0)
                prev_pos += 1
            else:
                out_list.append(line[4])
                prev_pos += 1
        return out_list
        
    def get_c(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        prev_pos = start 
        for line in _tabix_out:
            while prev_pos < int(line[1]) -1:
                out_list.append(0)
                prev_pos += 1
            else:
                out_list.append(line[5])
                prev_pos += 1
        return out_list
        
    def get_g(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        prev_pos = start 
        for line in _tabix_out:
            while prev_pos < int(line[1]) -1:
                out_list.append(0)
                prev_pos += 1
            else:
                out_list.append(line[6])
                prev_pos += 1
        return out_list
        
    def get_t(self,chr,start,end):
        _tabix_out = self.get_range(str(chr),start,end)
        out_list = []
        prev_pos = start 
        for line in _tabix_out:
            while prev_pos < int(line[1]) -1:
                out_list.append(0)
                prev_pos += 1
            else:
                out_list.append(line[7])
                prev_pos += 1
        return out_list
    """

# tools #

def sort_delim(delim_obj,col):
    s = delim.data[col]
    return sorted(range(len(s)),key=lambda i: s[i])

#def ordered_split(str,charlist):
#    out_list = [str]
#    for sc in charlist:
#        out_list = out_list.split(sc)
