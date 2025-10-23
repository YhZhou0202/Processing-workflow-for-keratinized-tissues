#!/bin/bash

### Trimming raw reads to generate clean reads
java -jar path/to/trimmomatic-0.39.jar PE -phred33 -threads 10 -trimlog logfilead raw.read1.gz raw.read2.gz clean.read1.gz out.trim1.gz clean.read2.gz out.trim2.gz ILLUMINACLIP:./MGI-adPE.fa:2:15:7:1:true SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3

###quality control
fastqc clean.read1.gz clean.read2.gz -o fq_quality

### Filter homologous sequences, generate comparison file unclassified.sam
bwa aln -l 20 -k 2 -t 2 "$ref" ./clean.read1.gz -f ./B.sa1.sai
bwa aln -l 20 -k 2 -t 2 "$ref" ./clean.read2.gz -f ./B.sa2.sai

bwa sampe -o 10000 -r '@RG\tID:B\tSM:B_SP\tPL:UNKNOWN\tLB:adapt' "$ref" ./B.sa1.sai ./B.sa2.sai ./clean.read1.gz ./clean.read2.gz -f ./unclassified.sam

###Modelling of DNA damage
damageprofiler -i unclassified.sam -o dam -r "$ref" -only_merged



###mtDNA assembly###
samtools fastq -N -F 4 -1 ./mapped1.fq -2 ./mapped2.fq unclassified.sam

gzip -f ./mapped1.fq

zcat ./mapped1.fq.gz | awk 'NR % 4 == 1 { sub(/\/1$/, "", $0); print }' > ./ID.txt
cat ./ID.txt | sed '/^@/ s/$/\/2/' > ./2ID.txt

zcat clean.read2.gz  | awk 'NR==FNR{a[$0]; next} $0 in a {flag=4} flag && flag--' ./2ID.txt - | gzip > ./mapped2.fq.gz

get_organelle_from_reads.py -1 mapped1.fq.gz -2 mapped2.fq.gz -s "$ref" -R 20 -t "$t" -F animal_mt --reduce-reads-for-coverage inf --max-reads inf -k 21,33,41 -o output