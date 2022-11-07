### Download data

wget ....

### Create file name list
```bash
data_directory=./data/
ls $data_directory/raw_data | cut -d '_' -f 1,2 |sort -u > ./sampleNames
```
### Quality trimming of samples
```bash
while read sampleNames; do (./scripts/quality_trimming.sh -s $sampleNames -i ./data/raw_data -o ./data/processed/ -p 30) done < sampleNames
```
### Get reference files
```bash
mkdir ./reference_files

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/685/GCF_011064685.1_ZJU_Sfru_1.0/GCF_011064685.1_ZJU_Sfru_1.0_genomic.fna.gz -P ./reference_files

gunzip ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/685/GCF_011064685.1_ZJU_Sfru_1.0/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gtf.gz -P ./reference_files

gunzip ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gtf.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/685/GCF_011064685.1_ZJU_Sfru_1.0/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gff.gz -P ./reference_files

gunzip ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/485/GCF_000838485.1_ViralProj14023/GCF_000838485.1_ViralProj14023_genomic.fna.gz -P ./reference_files

gunzip ./reference_files/GCF_000838485.1_ViralProj14023_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/485/GCF_000838485.1_ViralProj14023/GCF_000838485.1_ViralProj14023_genomic.gtf.gz -P ./reference_files

gunzip ./reference_files/GCF_000838485.1_ViralProj14023_genomic.gtf.gz

##insert how to get plasmid sequence

```
### Prepare a reference STAR index for rRNA removal
Due to the size of the data and a significant amount of rRNA that wasn't removed during library preparation this intermediate step is necessary to reduce data size. This step is not needed if rRNA was removed more successfully, or less data was acquired.

```bash
./scripts/prepare_star_rrna_index.sh -a ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gff -f ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.fna -p 30 -d ./reference_files
```

### Remove rRNA from Samples
This step shouldn't be run with more than 16 threads.
```bash
mkdir ./data/filtered_data
while read sampleNames; do (./scripts/run_rrna_star.sh -s $sampleNames -i ./data/processed/paired -p 16 -g ./reference_files/star_index_rrna -o ./data/filtered_data) done < sampleNames
```
### Create a hybrid genome reference
Using this code we create a hybrid reference sequence including references for Sf9, baculovirus and transgene sequences. 

```bash
cat ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.fna ./reference_files/GCF_000838485.1_ViralProj14023_genomic.fna ./reference_files/plasmid.fa >> ./overall.fasta

cat ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_genomic.gtf ./reference_files/GCF_000838485.1_ViralProj14023_genomic.gtf ./reference_files/plasmids.gtf >> ./overall.gtf
```

### Create a STAR index
```bash
./scripts/prepare_reference_star_index.sh -a ./reference_files/overall.gtf -f ./reference_files/overall.fasta -p 32 -d ./reference_files
```
### Map reads to Genome
```bash
mkdir ./data/mapped_data

while read sampleNames; do (./scripts/star_mapping.sh -s $sampleNames -i ./data/filtered_data/ -p 16 -g ./reference_files/star_index -o ./data/mapped_data) done < sampleNames
```

### Count mapped reads for differential expression analysis
```bash
while read sampleNames; do (./scripts/htseq_count.sh -s $sampleNames -g ./reference_files/overall.gtf -m ./data/mapped_data/ -o ./data/) done < sampleNames
```
### Separate counts
```bash
mkdir ./data/sf9_counts
mkdir ./data/bacv_counts

while read sampleNames; do(grep "ACNV" ./data/counts/"$sampleNames".count >> ./data/bacv_counts/"$sampleNames"_bacv.counts) done < sampleNames
while read sampleNames; do(grep "ACNV" -v  ./data/counts/"$sampleNames".count >> ./data/sf9_counts/"$sampleNames"_sf9.counts) done < sampleNames
```

### Run trinotate
This analysis is based on the instructions from the trinotate github page: https://github.com/Trinotate/Trinotate.github.io/wiki.
In order to download necessary programs, please refer to the tutorial.

#### Download additional files needed
```bash
mkdir ./trinotate
cd ./trinotate
mkdir ./reference_files
ref_dir=./reference_files

# Download protein reference data from RefSeq

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/685/GCF_011064685.1_ZJU_Sfru_1.0/GCF_011064685.1_ZJU_Sfru_1.0_cds_from_genomic.fna.gz -P $ref_dir

gunzip $ref_dir/GCF_011064685.1_ZJU_Sfru_1.0_cds_from_genomic.fna.gz 

# Prepare FASTA header for trinotate input
cat $ref_dir/GCF_011064685.1_ZJU_Sfru_1.0_cds_from_genomic.fna | sed 's/>.*protein_id=/>/' | sed 's/].*]//' >> $ref_dir/trinotate_sf9.fasta
```
#### Run trinotate to download necessary reference files for searches
```bash
./Trinotate-Trinotate-v3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

mv ./Pfam-A.hmm.gz $ref_dir
mv ./uniprot_sprot.dat.gz $ref_dir
mv ./uniprot_sprot.pep $ref_dir
mv ./Trinotate.sqlite $ref_dir

gunzip $ref_dir/Pfam-A.hmm.gz
```
#### Prepare databases for searches
```bash
# Prepare protein database for blast search 
makeblastdb -in $ref_dir/uniprot_sprot.pep -dbtype prot

# Prepare Pfam database 
hmmpress $ref_dir/Pfam-A.hmm

# Run Transdecoder 
./TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -m 10 -t $ref_dir/trinotate_sf9.fasta
```
#### Perform database searches
```bash
# Search reference transcripts
blastx -query $ref_dir/trinotate_sf9.fasta -db $ref_dir/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > $ref_dir/blastx.outfmt6

# Search transdecoder sequences
blastp -query ./trinotate_sf9.fasta.transdecoder_dir/longest_orfs.pep  -db $ref_dir/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3> $ref_dir/blastp.outfmt6

# Running HAMMER to identify protein domains
hmmscan --cpu 12 --domtblout TrinotatePFAM.out $ref_dir/Pfam-A.hmm ./trinotate_sf9.fasta.transdecoder_dir/longest_orfs.pep > pfam.log

# Run signalp-4.1 to predict signal peptides 
./signalp-4.1/signalp -f short -n $ref_dir/signalp.out ./trinotate_sf9.fasta.transdecoder_dir/longest_orfs.pep

./tmhmm-2.0c/bin/tmhmm --short < ./trinotate_sf9.fasta.transdecoder_dir/longest_orfs.pep > $ref_dir/tmhmm.out

./Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome $ref_dir/trinotate_sf9.fasta --path_to_rnammer ./RNAMMER/rnammer
```
#### Load genereted reports into Trinotate SQLite database
```bash
awk 'sub(/^>/, "")' ./reference_files/trinotate_sf9.fasta >> trinotate_headers.txt
cat ./reference_files/GCF_011064685.1_ZJU_Sfru_1.0_cds_from_genomic.fna | sed 's/>.*gene=/>/' | sed 's/].*]//' >> ./reference_files/trinotate_gene_sf9.fasta
awk 'sub(/^>/, "")' ./reference_files/trinotate_gene_sf9.fasta >> trinotate_gene_headers.txt
paste ./trinotate_gene_headers.txt ./trinotate_headers.txt | column -s '\t' -t >> ./trinotate_gene_trans_map.txt

./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite init --gene_trans_map ./trinotate_gene_trans_map.txt --transcript_fasta ./reference_files/trinotate_sf9.fasta --transdecoder_pep ./trinotate_sf9.fasta.transdecoder_dir/longest_orfs.pep

./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite LOAD_swissprot_blastp ./reference_files/blastp.outfmt6 
./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite LOAD_swissprot_blastx ./reference_files/blastx.outfmt6 
./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite LOAD_pfam ./reference_files/TrinotatePFAM.out
./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite LOAD_tmhmm ./tmhmm.out
./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite LOAD_signalp ./signalp.out

mkdir ./trinotate_output

./Trinotate-Trinotate-v3.2.2/Trinotate ./reference_files/Trinotate.sqlite report > ./trinotate_output/trinotate_sf9_annotation_report.xls


./Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ./trinotate_output/trinotate_sf9_annotation_report.xls -G --include_ancestral_terms > ./trinotate_output/sf9_go_annotations.txt


./Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ./trinotate_output/trinotate_sf9_annotation_report.xls -T --include_ancestral_terms > ./trinotate_output/sf9_go_annotations_v2.txt

./Trinotate-Trinotate-v3.2.2/util/Trinotate_get_feature_name_encoding_attributes.pl ./trinotate_output/trinotate_sf9_annotation_report.xls > ./trinotate_output/annot_feature_map.txt

```
