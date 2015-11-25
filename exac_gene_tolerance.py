#!/usr/bin/env python2.7

genes =[]
vc1 = ["ESS","FS","SG"]
vc2 = ["NSY","SS5","IF","IM","SL","EE"]
vc3 = ["SY","SS","INT","5PU","3PU"]
vclass1 = {}
vclass2 = {}
vclass3 = {}

with open("/scratch/cancgene/gst/databases/20150926_ExAC.r0.3.sites_cavaOutput_uniqENSTandCSN.txt","r") as FH:
    head= FH.readline()
    for line in FH:
        l=line.rstrip()
        exac_var = l.split("\t")
        if exac_var[4] in genes:
            if exac_var[8] in vc1:
                vclass1[exac_var[4]] += int(exac_var[12])
            elif exac_var[8] in vc2:
                vclass2[exac_var[4]] += int(exac_var[12])
            elif exac_var[8] in vc3:
                vclass3[exac_var[4]] += int(exac_var[12])
        else:
            if exac_var[4] == ".":  
                continue
            genes.append(exac_var[4])
            if exac_var[8] in vc1:
                vclass1[exac_var[4]] = int(exac_var[12])
                vclass2[exac_var[4]] = 0
                vclass3[exac_var[4]] = 0
            elif exac_var[8] in vc2:
                vclass1[exac_var[4]] = 0
                vclass2[exac_var[4]] = int(exac_var[12])
                vclass3[exac_var[4]] = 0
            elif exac_var[8] in vc3:
                vclass1[exac_var[4]] = 0
                vclass2[exac_var[4]] = 0
                vclass3[exac_var[4]] = int(exac_var[12])
            else:
                vclass1[exac_var[4]] = 0
                vclass2[exac_var[4]] = 0
                vclass3[exac_var[4]] = 0
out = open("exac_gene_variant_counts.txt","w")
for gene in set(genes):
    outline= "\t".join([gene,str(vclass1[gene]),str(vclass2[gene]),str(vclass3[gene])])+"\n"
    out.write(outline)