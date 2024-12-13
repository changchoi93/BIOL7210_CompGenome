#!/usr/bin/env sh

#Fetch the genome mentioned in the instructions and unzip
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/427/945/GCF_034427945.1_ASM3442794v1/GCF_034427945.1_ASM3442794v1_genomic.fna.gz
pigz -d GCF_034427945.1_ASM3442794v1_genomic.fna.gz

#I chose Prodigal as my package to use so activate an environment with the package installed
conda activate ex4_pt2

#Run Prodigal, storing stderr and stdout as a single logfile
prodigal -i ~/ex4/GCF_034427945.1_ASM3442794v1_genomic.fna -c -m -f gff -o cds.gff 2>&1 | tee log.txt

#I chose barrnap as my package to use so activate an environment with the package installed
conda deactivate
conda activate ex4_pt1

#Run barrnap and extract all 16S rRNA sequences
barrnap GCF_034427945.1_ASM3442794v1_genomic.fna | grep "Name=16S_rRNA;product=16S ribosomal RNA" > 16S.gff
bedtools getfasta -fi GCF_034427945.1_ASM3442794v1_genomic.fna -bed 16S.gff -fo 16S.fa
conda deactivate
conda activate ex4_pt2
pigz -9f 16S.fa cds.gff log.txt

#Use the NCBI Blast website to identify the top 5 hits, use blastn because we are querying nucleotides, use the rRNA/ITS database (also select 16S Bacteria), and select the top 5 hits by max score
#Use the ompA.faa provided in the repo to predict genes and do research to write functional descriptions of each protein in a Word document