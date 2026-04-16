#!/bin/bash

INPUT=$1
WINDOW=$2
BASE="$(basename $INPUT | sed 's/\_rmdup.bam$//g')"
OUTPUT=$3

SCALE=$(echo "1000000/$(samtools view -c $INPUT)" | bc -l)

bedtools intersect -c -a $WINDOW -b $INPUT | \
awk -v scale=$SCALE '{print $1,$2,$3,$4*1000*(scale/($3-$2)) }' OFS='\t' > $OUTPUT

echo $OUTPUT
