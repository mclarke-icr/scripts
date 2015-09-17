#!/bin/bash
DATE=$( date +%Y%m%d_%H.%M.%S )
cd /scratch/cancgene/mclarke/
git add *
git add lib/*
git commit -m "autocommit $DATE"
git push origin master
