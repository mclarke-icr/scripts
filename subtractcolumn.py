#!/usr/bin/env python

import sys

first = [ line.rstrip() for line in open(sys.argv[1],"r") ]
fc = int(sys.argv[2])
sc = int(sys.argv[4])
second = []
for line in open(sys.argv[3],"r"):
    l= line.rstrip()
    cols = l.split("\t")
    second.append(cols[sc])


for f in first:
    f_cols = f.split("\t")
    if f_cols[fc] in second:
        pass
    else:
        print f