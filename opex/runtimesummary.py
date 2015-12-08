#!/usr/bin/env python2.7

import sys
from datetime import datetime

filelist = sys.argv[1]

files = [ line.rstrip() for line in open(filelist,"r") ]
intervals ={}
for f in files:
    logf = open(f,"r")
    sample = "_".join(f.split("_")[:-2])
    timepoints = []
    for line in logf:
        fields = line.split(" ")
        (day,month,year) = fields[0].split(".")
        (hour,minute,second) = fields[1].split(":")
        dt = datetime(int(year),int(month),int(day),int(hour),int(minute),int(second))
        timepoints.append(dt)
        #print dt.isoformat()
    timediff = []
    for t_idx in range(len(timepoints)):
        if t_idx == 0:
            continue
        t_delta = timepoints[t_idx] - timepoints[t_idx-1]
        timediff.append(t_delta)
    intervals[sample] = timediff
    logf.close()

count = 0
sum = {}
out = open("all_opex_log.txt","w")
for samp in intervals:
    count +=1
    out_times = []
    for t_delta in intervals[samp]:
        out_times.append(str(int(t_delta.total_seconds())))
    out_line=  samp+"\t"+"\t".join(out_times)+"\n"
    out.write(out_line)
