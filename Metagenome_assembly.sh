#!/bin/bash
####Metagenome Assembly####

#Part 1 - Environment set-up

#create and activate Conda environment
localdisk/software/anaconda3/bin/conda env create -n NGG4 --file /localdisk/data/NGG/conda_envs/NGG4.yml
source /localdisk/software/anaconda3/bin/activate NGG4

#link data
ln -s /localdisk/software/blobtools-1.1/data/nodesDB.txt ./nodesDB.txt
ln -s /localdisk/data/NGG/ERR260505_1.fastq.gz ./ERR260505_1.fastq.gz
ln -s /localdisk/data/NGG/ERR260505_2.fastq.gz ./ERR260505_2.fastq.gz
ln -s /localdisk/data/NGG/NexteraPE-PE.fa ./NexteraPE-PE.fa
ln -s /localdisk/data/NGG/taxid_map ./taxid_map
ln -s /localdisk/data/NGG/genomes.fasta ./g

#Part 2 - Quality control and read trimming

#Quality control the raw data using FastQC
fastqc -t 2 ERR260505_1.fastq.gz ERR260505_2.fastq.gz
firefox ERR260505_1_fastqc.html ERR260505_2_fastqc.html &

#Trim and filter raw reads using Skewer
mkdir trimmed_reads
skewer-0.2.2-linux-x86_64 -n -Q 20 -l 75 -t 2 -m any -x NexteraPE-PE.fa ERR260505_1.fastq.gz ERR260505_2.fastq.gz -o trimmed_reads/ERR260505

#Part 3 - Assembly using SPAdes

#Assembly using SPAdes in standard genome mode
spades.py --only-assembler -m 30 -t 2 -1 trimmed_reads/ERR260505-trimmed-pair1.fastq -2 trimmed_reads/ERR260505-trimmed-pair2.fastq -o spades

#examine the SPAdes assembly including N50 and othe assembly metrics:
~/scaffold_stats.pl -f spades/scaffolds.fasta -t 100 -h

#Assembly using SPAdes in metagenome mode
spades.py --meta --only-assembler -m 30 -t 4 -1 trimmed_reads/ERR260505-trimmed-pair1.fastq -2 trimmed_reads/ERR260505-trimmed-pair2.fastq -o spades_meta
~/scaffold_stats.pl -f spades_meta/scaffolds.fasta -t 100 -h

#Part 4 - Compare metagenome assemblies using MetaQUAST
mv spades/scaffolds.fasta spades/spades.fasta
mv spades_meta/scaffolds.fasta spades_meta/spades_meta.fasta

metaquast -o metaquast/ --no-plots -t 4 spades/spades.fasta spades_meta/spades_meta.fasta
firefox metaquast/report.html metaquast/krona_charts/summary_taxonomy_chart.html &

#Part 5 - Visualise the taxonomic content of metagenome assemblies using blobtools

#BLAST annotation
makeblastdb -taxid_map taxid_map -parse_seqids -in genomes.fasta -dbtype nucl
blastn -query spades_meta/spades_meta.fasta -db genomes.fasta -outfmt '6 qseqid staxids bitscore std' -num_threads 2 > blast.out

#Align raw reads to the metagenome assembly using BWA

#Create an index of the SPAdes assembly generated using the metagenome mode:
bwa index -a is spades_meta/spades_meta.fasta

#Find the suffix array coordinates of the single-end reads:
mkdir aligned_reads
bwa aln -t 4 spades_meta/spades_meta.fasta trimmed_reads/ERR260505-trimmed-pair1.fastq > aligned_reads/ERR260505_1.sai
bwa aln -t 4 spades_meta/spades_meta.fasta trimmed_reads/ERR260505-trimmed-pair2.fastq > aligned_reads/ERR260505_2.sai

#Generate paired-end alignments in SAM format:
bwa sampe spades_meta/spades_meta.fasta aligned_reads/ERR260505_1.sai aligned_reads/ERR260505_2.sai trimmed_reads/ERR260505-trimmed-pair1.fastq trimmed_reads/ERR260505-trimmed-pair2.fastq > aligned_reads/ERR260505.sam
samtools sort aligned_reads/ERR260505.sam -o aligned_reads/ERR260505.sorted.bam
samtools index aligned_reads/ERR260505.sorted.bam

#Draw blob plots using blobtools
blobtools create -i spades_meta/spades_meta.fasta -b aligned_reads/ERR260505.sorted.bam -t blast.out --db nodesDB.txt
blobtools plot -i blobDB.json -m -r species
mkdir blobplots
mv *.png blobplots/