#!/usr/bin/env bash
#### Count reads aligned to reference genome for differential expression analysis
#### Inputs are: 1) sample names, 2) bam input directory, 3) reference GTF and 4) output directory
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample name"
        echo "-m = Mapping directory"
        echo "-g = Reference GTF"
        echo "-o = path for output count files"
        exit 2
fi
while getopts s:m:g:o: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      g) REF_GTF=${OPTARG};;
      m) MAP_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

mkdir -p $OUT_DIR/counts
echo $SAMPLE_ID
htseq-count -n 10 -r pos -f bam -t exon -i gene_id -s reverse $MAP_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam $REF_GTF > "$OUT_DIR"/counts/"$SAMPLE_ID".count

# END
