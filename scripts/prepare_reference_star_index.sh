#!/usr/bin/env bash
#### Make reference index for STAR
#### inputs are: 1) Hybrid GTF and 2) Hybrid FASTA 3) NUM processors and 4) output genome directory
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

while getopts a:f:p:d: option
  do
    case "${option}"
      in
      a) GTF=${OPTARG};;
      f) FASTA=${OPTARG};;
      p) THREADS=${OPTARG};;
      d) GENOME_DIR=${OPTARG};;
    esac
done

mkdir $GENOME_DIR/star_index
#Create a STAR index

STAR --runMode genomeGenerate \
     --sjdbOverhang 124 \
     --genomeChrBinNbits 16 \
     --genomeDir $GENOME_DIR/star_index \
     --genomeFastaFiles $FASTA \
     --genomeSAindexNbases 13 \
     --sjdbGTFfile $GTF \
     --runThreadN $THREADS

# END
