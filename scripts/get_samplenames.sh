#!/usr/bin/env bash
### Get sequencing codes and statistics for further evaluation
### Inputs are: i), o)
### Written by iBET: nikolaus.virgolini@ibet.pt

if(($# == 0)); then
	echo "Usage:"
	echo "-i = input directory"
	echo "-o = output_directory"
	exit 2
fi

while getopts i:o option
	do
	  case "${option}"
	     in
	     i) IN_DIR=${OPTARG};;
	     o) OUT_DIR=${OPTARG};;
	  esac
done

###Make a file name list
ls $IN_DIR/raw_data | sort -u > ./sampleNames
###Count raw reads per file
while read sampleNames; do \
    current_forward_reads=$(ls $IN_DIR/ | grep $sampleNames | head  -n 1); 
    current_reverse_reads=$(ls $IN_DIR/ | grep $sampleNames | tail  -n 1); 
    num_forward_reads=$(zcat "$$IN_DIR/$current_forward_reads" | echo $((`wc -l`/4)));
    num_reverse_reads=$(zcat "$IN_DIR/$current_reverse_reads" | echo $((`wc -l`/4)));
    echo -e $sampleName'\t'$num_forward_reads'\t'$num_reverse_reads >> $OUT_DIR/raw.data.counts \ 
done < ./sampleNames

#END
