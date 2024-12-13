#!/usr/bin/env bash

conda activate ex3

mkdir ~/ex6/Raw_FastQs
cd ~/ex6/Raw_FastQs
for accession in SRR1993270 SRR1993271 SRR1993272 SRR2984947 SRR2985018 SRR3214715 SRR3215024 SRR3215107 SRR3215123 SRR3215124; do
  fasterq-dump \
   "${accession}" 
   --outdir . \
   --split-files \
   --skip-technical
done
pigz -9 *.fastq

fastp -i /home/chris/ex6/Raw_FastQs/SRR1993270_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR1993270_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR1993270.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR1993270.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR1993270.json --html ~/ex6/Cleaned_FastQs/SRR1993270.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR1993271_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR1993271_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR1993271.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR1993271.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR1993271.json --html ~/ex6/Cleaned_FastQs/SRR1993271.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR1993272_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR1993272_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR1993272.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR1993272.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR1993272.json --html ~/ex6/Cleaned_FastQs/SRR1993272.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR2984947_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR2984947_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR2984947.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR2984947.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR2984947.json --html ~/ex6/Cleaned_FastQs/SRR2984947.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR2985018_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR2985018_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR2985018.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR2985018.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR2985018.json --html ~/ex6/Cleaned_FastQs/SRR2985018.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR3214715_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR3214715_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR3214715.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR3214715.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR3214715.json --html ~/ex6/Cleaned_FastQs/SRR3214715.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR3215024_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR3215024_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR3215024.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR3215024.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR3215024.json --html ~/ex6/Cleaned_FastQs/SRR3215024.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR3215107_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR3215107_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR3215107.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR3215107.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR3215107.json --html ~/ex6/Cleaned_FastQs/SRR3215107.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR3215123_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR3215123_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR3215123.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR3215123.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR3215123.json --html ~/ex6/Cleaned_FastQs/SRR3215123.html
fastp -i /home/chris/ex6/Raw_FastQs/SRR3215124_1.fastq.gz -I /home/chris/ex6/Raw_FastQs/SRR3215124_2.fastq.gz -o ~/ex6/Cleaned_FastQs/SRR3215124.R1.fq.gz -O ~/ex6/Cleaned_FastQs/SRR3215124.R2.fq.gz --json ~/ex6/Cleaned_FastQs/SRR3215124.json --html ~/ex6/Cleaned_FastQs/SRR3215124.html

mkdir ~/ex6/Assemblies
for read in ~/ex6/Cleaned_FastQs/*.R1.fq.gz; do
  sample="$(basename ${read} .R1.fq.gz)"
  skesa \
   --reads "${read}","${read%R1.fq.gz}R2.fq.gz" \
   --cores 4 \
   --min_contig 1000 \
   --contigs_out ~/ex6/Assemblies/"${sample}".fna
done

cd ~/ex6/Assemblies
ls -alhS *.fna

conda deactivate

conda create -n harvestsuite -c bioconda parsnp harvesttools figtree -y
conda activate harvestsuite

mkdir ~/ex6/parsnp_input_assemblies
cd ~/ex6/parsnp_input_assemblies

for file in ~/ex6/Assemblies/*.fna; do
    cp "$file" ~/ex6/parsnp_input_assemblies
done

ls -alhtr *.{fa,fna}

cd ~/ex6

parsnp \
 -d parsnp_input_assemblies \
 -r ~/ex6/parsnp_input_assemblies/SRR3214715.fna \
 -o parsnp_outdir \
 -p 4

figtree parsnp_outdir/parsnp.tree

#From the pop-up GUI, save the tree as a PDF (or another format of your choice).

mv parsnp.tree.pdf Choi.pdf

#Look at the tree and determine which SRA accessions were outliers.

echo 'SRR1993272 SRR1993271 SRR1993270 SRR2984947 SRR2985018' > Choi.txt