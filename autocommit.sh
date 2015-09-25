#!/bin/bash
COMMENT=$1
DATE=$( date +%Y%m%d_%H.%M.%S )
cd /scratch/cancgene/mclarke/scripts
git add *
git add lib/*
git add opex/*
if [ -z "$COMMENT" ]
then
  git commit -m "$DATE autocommit"
else
  git commit -m "$DATE $COMMENT"
fi
git push origin master
