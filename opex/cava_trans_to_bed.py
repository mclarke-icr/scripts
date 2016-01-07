#!/bin/env python

exome_ref = open('/scratch/cancgene/mclarke/opex/opex-v1.0.0_old/exome_65_GRCh37','r')
out=open('/scratch/cancgene/mclarke/opex/exome_65_GRCh37.bed','w')

for l in exome_ref:
    line = l.rstrip()
    fields = line.split("\t")
    boundries = []
    strand='0'
    if fields[3] == '1':
        boundries = fields[7:]
        strand='+'
    elif fields[3] == '-1':
        boundries = fields[7:]
        strand='-'
    else:
        print "nostrand"

    for x in range(0,len(boundries),2):
        bedout = "\t".join([fields[2],boundries[x],boundries[x+1],fields[0],"0",strand])+"\n"
        out.write(bedout)
        