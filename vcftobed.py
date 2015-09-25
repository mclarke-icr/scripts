#!/bin/env python

import argparse



###

parser = argparse.ArgumentParser(description='gives vcftobed -i infile.vcf -o outfile')
parser.add_argument("-i","--input",required=True)
parser.add_argument("-f","--family")
parser.add_argument("-o","--output")
args = vars(parser.parse_args())
infile = args["input"]
family = args["family"]
outfile = args["output"]