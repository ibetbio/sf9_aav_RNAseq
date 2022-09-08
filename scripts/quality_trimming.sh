#!/usr/bin/env bash
#### Trim low quality bases using trimmomatic
#### inputs are: 1) sample ID, 2) the data directory, 3) the  output directory and 4) the number of threads
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-o = output directory"
        echo "-p = num threads"
        exit 2
fi

while getopts s:i:o:p: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      i) IN_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
      p) NUM_THREADS=${OPTARG};;
    esac
done

mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/paired $OUT_DIR/unpaired
trim_directory=/mnt/HDD2/nikolaus/Trimmomatic-0.39/

java -jar  $trim_directory/trimmomatic-0.39.jar PE \
     -threads $NUM_THREADS \
     $IN_DIR/"$SAMPLE_ID"_R1.fastq.gz $IN_DIR/"$SAMPLE_ID"_R2.fastq.gz \
     $OUT_DIR/paired/"$SAMPLE_ID"_R1.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_R1.fastq.gz \
     $OUT_DIR/paired/"$SAMPLE_ID"_R2.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_R2.fastq.gz \
     SLIDINGWINDOW:4:20 \
     MINLEN:36 \
     -trimlog $OUT_DIR/"$SAMPLE_ID".trimmomatic.log
# END
