#!/usr/bin/env sh

#activate my previously created conda environment with entrez-direct, sra-tools, fastqc, trimmomatic, skesa, spades, and pigz installed
conda activate ex3

#make a raw data folder, fetch the SRR26244988 data using fasterq-dump, then compress it
mkdir -pv ~/biol7210/ex3/raw_data
fasterq-dump SRR26244988 --threads 1 --outdir ~/biol7210/ex3/raw_data --split-files --skip-technical
pigz -9f ~/biol7210/ex3/raw_data/*.fastq

#look at the read quality using fastqc
mkdir -v ~/biol7210/ex3/raw_qa
fastqc --threads 2 --outdir ~/biol7210/ex3/raw_qa ~/biol7210/ex3/raw_data/SRR26244988_1.fastq.gz ~/biol7210/ex3/raw_data/SRR26244988_2.fastq.gz

#looking at the fastqc report, the read qualities were mostly good so trim the reads using the same parameters as the example given in class (dropping the read if the average quality is below 30 and using a sliding window trimming approach)
mkdir -v ~/biol7210/ex3/trim
cd ~/biol7210/ex3/trim
trimmomatic PE -phred33 ~/biol7210/ex3/raw_data/SRR26244988_1.fastq.gz ~/biol7210/ex3/raw_data/SRR26244988_2.fastq.gz ~/biol7210/ex3/trim/r1.paired.fq.gz ~/biol7210/ex3/trim/r1_unpaired.fq.gz ~/biol7210/ex3/trim/r2.paired.fq.gz ~/biol7210/ex3/trim/r2_unpaired.fq.gz SLIDINGWINDOW:5:30 AVGQUAL:30 1> trimmo.stdout.log 2> trimmo.stderr.log
cat ~/biol7210/ex3/trim/r1_unpaired.fq.gz ~/biol7210/ex3/trim/r2_unpaired.fq.gz > ~/biol7210/ex3/trim/singletons.fq.gz
rm -v ~/biol7210/ex3/trim/*unpaired*

#do genome assembly with spades
mkdir -v ~/biol7210/ex3/asm
cd ~/biol7210/ex3/asm
spades.py -1 ~/biol7210/ex3/trim/r1.paired.fq.gz -2 ~/biol7210/ex3/trim/r2.paired.fq.gz -o ~/biol7210/ex3/asm

#switch to an environment with python 2.7 and biopython
conda deactivate
conda activate bpy2

#use the filter.contigs.py script (which was provided) with the double the default parameters for minimum coverage and length
~/biol7210/ex3/filter.contigs.py -i ~/biol7210/ex3/asm/contigs.fasta -o ~/biol7210/ex3/asm/filtered_assembly.fna -c 10 -l 1000