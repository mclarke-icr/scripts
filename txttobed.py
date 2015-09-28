#!/bin/env python

import argparse

# inputs
parser = argparse.ArgumentParser(description='gives vcftobed -i infile.txt -o outfile')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-o","--output")
args = vars(parser.parse_args())
infile = args["input"]
outfile = args["output"]
#outfile
if outfile == None:
    outfile= "".join(infile.split(".txt")) + ".bed"
out=open(outfile,"w")
#loop

with open(infile, 'r') as f:
    for line in  f.readlines():
        if line.startswith('#'):
            continue
        data = line.split("\t")
        out.write("\t".join([data[0],str(int(data[1])-1),data[1],data[7]])+"\n")