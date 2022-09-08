#!/usr/bin/env bash
#### Map reads to reference rRNA index to remove rRNA
#### The inputs are: 1) sample ID and 2) the data directories 3) the output directory, and  4) number of threads
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-p = num processors"
        echo "-g = path to star genome index"
        echo "-o = path for output files"
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
--genomeLoad NoSharedMemory \
--seedSearchStartLmaxOverLread .5 \
--outFilterMultimapNmax 1000 \
--outFilterMismatchNmax 2 \
--genomeDir $INDEX \
--outFileNamePrefix $OUTDIR/"$SAMPLEID" \
--readFilesIn $INDIR/"$SAMPLEID"_R1.fastq.gz $INDIR/"$SAMPLEID"_R2.fastq.gz \
--readFilesCommand zcat \
--outReadsUnmapped Fastx

gzip $OUTDIR/"$SAMPLEID"Unmapped.out.mate1
gzip $OUTDIR/"$SAMPLEID"Unmapped.out.mate2

#In order to save disc space we delete the SAM file after mapping. This step is optional
rm $OUTDIR/"$SAMPLEID"Aligned.out.sam

# END
