#!/usr/bin/env python

import sys

first = [ line.rstrip() for line in open(sys.argv[1],"r") ]
second = [ line.rstrip() for line in open(sys.argv[2],"r") ]

for f in first:
    if f in second:
        pass
    else:
        print f