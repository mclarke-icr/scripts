#!/bin/bash
SS=COG1670_trio.txt
module load python/2.7.5
cd /scratch/cancgene/mclarke/opex/COG1670_fromVCF
# selecting varaints unseen in parents
../../scripts/opex/opex_denovo.py -i COG1670_trio.txt

#running coverview
SAMP=( $( cat $SS | awk '{print $1}' ) )
FAM=( $( cat $SS | awk '{print $2}' ) )
REL=( $( cat $SS | awk '{print $3}' ) )
END=( $( echo $( cat $SS | wc -l ) -1 | bc ) )
for i in $( seq 1 $END) # starts at 1 due to header
do
  python2.7 ../opex-v1.0.0/tools/CoverView-v1.1.0/CoverView.py -i ${SAMP[$i]}_picard.bam -b ${FAM[$i]}.bed -o ${SAMP[$i]}_coverview
done

# addibg covrage stats back to the .txt
python2.7 ../../scripts/opex/add_cov.py -i COG1670_trio_denovo.txt -f COG1670_trio.txt
