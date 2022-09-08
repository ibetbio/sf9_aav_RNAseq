#!/usr/bin/env bash
#### Map reads to reference file
#### The inputs are: 1) sample ID and 2) the data directories 3) the output directory, and  4) number of threads
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-p = num processors"
        echo "-g = path to star genome index"
        echo "-o = path for output BAMs"
        exit 2
fi
while getopts s:i:p:g:o: option
  do
    case "${option}"
      in
      s) SAMPLEID=${OPTARG};;
      i) INDIR=${OPTARG};;
      p) THREADS=${OPTARG};;
      g) INDEX=${OPTARG};;
      o) OUTDIR=${OPTARG};;
    esac
done

STAR \
--runThreadN $THREADS \
--genomeDir $INDEX \
--readFilesCommand gunzip -c \
--outFileNamePrefix $OUTDIR/"$SAMPLEID" \
--readFilesIn $INDIR/"$SAMPLEID"Unmapped.out.mate1.gz $INDIR/"$SAMPLEID"Unmapped.out.mate2.gz \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic

# END

