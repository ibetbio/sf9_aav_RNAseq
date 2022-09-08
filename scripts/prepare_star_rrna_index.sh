#!/usr/bin/env bash
#### Prepare reference file for rRNA removal and build STAR index
#### inputs are: 1) RefSeq GFF and 2) RefSeq FASTA 3) NUM processors and 4) reference genome directory
#### Written by iBET: nikolaus.virgolini@ibet.pt 2022-01

while getopts a:f:p:d: option
  do
    case "${option}"
      in
      a) GFF=${OPTARG};;
      f) FASTA=${OPTARG};;
      p) THREADS=${OPTARG};;
      d) REF_DIR=${OPTARG};;
    esac
done

mkdir $REF_DIR/rrna
mkdir $REF_DIR/star_index_rrna

#Grep all rRNA sequences from the GFF file
grep -i "rRNA" $GFF >> $REF_DIR/rrna/rrna.gff

#Create a reference fasta file for reating the STAR index
gffread -w $REF_DIR/rrna/transcripts.fa -g $FASTA $REF_DIR/rrna/rrna.gff

#Create a STAR index

STAR --runMode genomeGenerate \
     --genomeSAindexNbases 6 \
     --genomeDir $REF_DIR/star_index_rrna \
     --genomeFastaFiles $REF_DIR/rrna/transcripts.fa \
     --runThreadN $THREADS

# END
