#!/usr/bin/env sh

#Activate a conda environment with sra-tools, trimmomatic, skesa
conda activate ex3

#Download the 3 Campylobacter jejuni isolates using their SRA numbers
fasterq-dump SRR3214715 SRR3215024 SRR3215107 --threads 1 --outdir ~/hw3_gtqa/raw_data --split-files --skip-technical

#Read clean with trimmomatic
trimmomatic PE -phred33 ~/hw3_gtqa/raw_data/SRR3214715_1.fastq ~/hw3_gtqa/raw_data/SRR3214715_2.fastq ~/hw3_gtqa/trim/SRR3214715.r1.paired.fq ~/hw3_gtqa/trim/SRR3214715.r1.unpaired.fq ~/hw3_gtqa/trim/SRR3214715.r2.paired.fq ~/hw3_gtqa/trim/SRR3214715.r2.unpaired.fq SLIDINGWINDOW:5:30 AVGQUAL:30 1> trimmo.stdout.log 2> trimmo.stderr.log
trimmomatic PE -phred33 ~/hw3_gtqa/raw_data/SRR3215024_1.fastq ~/hw3_gtqa/raw_data/SRR3215024_2.fastq ~/hw3_gtqa/trim/SRR3215024.r1.paired.fq ~/hw3_gtqa/trim/SRR3215024.r1.unpaired.fq ~/hw3_gtqa/trim/SRR3215024.r2.paired.fq ~/hw3_gtqa/trim/SRR3215024.r2.unpaired.fq SLIDINGWINDOW:5:30 AVGQUAL:30 1> trimmo.stdout.log 2> trimmo.stderr.log
trimmomatic PE -phred33 ~/hw3_gtqa/raw_data/SRR3215107_1.fastq ~/hw3_gtqa/raw_data/SRR3215107_2.fastq ~/hw3_gtqa/trim/SRR3215107.r1.paired.fq ~/hw3_gtqa/trim/SRR3215107.r1.unpaired.fq ~/hw3_gtqa/trim/SRR3215107.r2.paired.fq ~/hw3_gtqa/trim/SRR3215107.r2.unpaired.fq SLIDINGWINDOW:5:30 AVGQUAL:30 1> trimmo.stdout.log 2> trimmo.stderr.log
rm -v ~/hw3_gtqa/trim/*unpaired*

#Assemble with skesa
skesa --reads ~/hw3_gtqa/trim/SRR3214715.r1.paired.fq ~/hw3_gtqa/trim/SRR3214715.r2.paired.fq --contigs_out ~/hw3_gtqa/asm/SRR3214715_skesa.fna 1> skesa.stdout.txt 2> skesa.stderr.txt
skesa --reads ~/hw3_gtqa/trim/SRR3215024.r1.paired.fq ~/hw3_gtqa/trim/SRR3215024.r2.paired.fq --contigs_out ~/hw3_gtqa/asm/SRR3215024_skesa.fna 1> skesa.stdout.txt 2> skesa.stderr.txt
skesa --reads ~/hw3_gtqa/trim/SRR3215107.r1.paired.fq ~/hw3_gtqa/trim/SRR3215107.r2.paired.fq --contigs_out ~/hw3_gtqa/asm/SRR3215107_skesa.fna 1> skesa.stdout.txt 2> skesa.stderr.txt

#Switch to an environment with python 2.7 and biopython
conda deactivate
conda activate bpy2

#Use Dr. Gulvik's filter.contigs.py to filter out low coverage and short contigs
~/biol7210/ex3/filter.contigs.py -i ~/hw3_gtqa/asm/SRR3214715_skesa.fna -o ~/hw3_gtqa/asm/SRR3214715_filtered.fna -c 10 -l 1000
~/biol7210/ex3/filter.contigs.py -i ~/hw3_gtqa/asm/SRR3215024_skesa.fna -o ~/hw3_gtqa/asm/SRR3215024_filtered.fna -c 10 -l 1000
~/biol7210/ex3/filter.contigs.py -i ~/hw3_gtqa/asm/SRR3215107_skesa.fna -o ~/hw3_gtqa/asm/SRR3215107_filtered.fna -c 10 -l 1000

#Download the reference type strain for Campylobacter jejuni to compare to, and then decompress the file
cd ~/hw3_gtqa/fastani
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/GCF_000009085.1_ASM908v1_genomic.fna.gz
gunzip -kv *.fna.gz

#Switch to an environment with fastani
conda deactivate
conda activate fastani

#Rename files for readability
cp -a ~/hw3_gtqa/asm/SRR3214715_filtered.fna ~/hw3_gtqa/fastani/
cp -a ~/hw3_gtqa/asm/SRR3215024_filtered.fna ~/hw3_gtqa/fastani/
cp -a ~/hw3_gtqa/asm/SRR3215107_filtered.fna ~/hw3_gtqa/fastani/
mv -v GCF_000009085.1_ASM908v1_genomic.fna GCF_000009085.1_ASM908v1.fna
mv -v SRR3214715_filtered.fna SRR3214715.fna
mv -v SRR3215024_filtered.fna SRR3215024.fna
mv -v SRR3215107_filtered.fna SRR3215107.fna

#Run fastani to compare similarity
fastANI --query SRR3214715.fna --ref GCF_000009085.1_ASM908v1.fna --output SRR3214715_fastani.tsv
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' SRR3214715_fastani.tsv > SRR3214715_fastani_alignment.tsv
sed "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" SRR3214715_fastani_alignment.tsv > SRR3214715_fastani_alignment_header.tsv

fastANI --query SRR3215024.fna --ref GCF_000009085.1_ASM908v1.fna --output SRR3215024_fastani.tsv
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' SRR3215024_fastani.tsv > SRR3215024_fastani_alignment.tsv
sed "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" SRR3215024_fastani_alignment.tsv > SRR3215024_fastani_alignment_header.tsv

fastANI --query SRR3215107.fna --ref GCF_000009085.1_ASM908v1.fna --output SRR3215107_fastani.tsv
awk '{alignment_percent = $4/$5*100} {alignment_length = $4*3000} {print $0 "\t" alignment_percent "\t" alignment_length}' SRR3215107_fastani.tsv > SRR3215107_fastani_alignment.tsv
sed "1i Query\tReference\t%ANI\tNum_Fragments_Mapped\tTotal_Query_Fragments\t%Query_Aligned\tBasepairs_Query_Aligned" SRR3215107_fastani_alignment.tsv > SRR3215107_fastani_alignment_header.tsv

#Rename the first two columns to just include the accession numbers as stated in the directions
awk 'BEGIN{FS=OFS="\t"} NR==1 {print; next} {for (i=1; i<=2; i++) sub(/\.fna$/, "", $i)} 1' SRR3214715_fastani_alignment_header.tsv > SRR3214715_fastani_final.tsv
awk 'BEGIN{FS=OFS="\t"} NR==1 {print; next} {for (i=1; i<=2; i++) sub(/\.fna$/, "", $i)} 1' SRR3215024_fastani_alignment_header.tsv > SRR3215024_fastani_final.tsv
awk 'BEGIN{FS=OFS="\t"} NR==1 {print; next} {for (i=1; i<=2; i++) sub(/\.fna$/, "", $i)} 1' SRR3215107_fastani_alignment_header.tsv > SRR3215107_fastani_final.tsv

#Combine the 3 tsv files into a single one for submission
(cat SRR3214715_fastani_final.tsv && awk 'NR > 1' SRR3215024_fastani_final.tsv && awk 'NR > 1' SRR3215107_fastani_final.tsv) > fastani.tsv

#Switch working directory and environment to one with mlst
conda deactivate
conda activate mlst
cd ~/hw3_gtqa/mlst

#Copy files over for ease (they're small enough)
cp -a ~/hw3_gtqa/fastani/SRR3214715.fna ~/hw3_gtqa/mlst/
cp -a ~/hw3_gtqa/fastani/SRR3215024.fna ~/hw3_gtqa/mlst/
cp -a ~/hw3_gtqa/fastani/SRR3215107.fna ~/hw3_gtqa/mlst/
cp -a ~/hw3_gtqa/fastani/GCF_000009085.1_ASM908v1.fna ~/hw3_gtqa/mlst/

#Run mlst to genotype and annotate
mlst SRR3214715.fna > SRR3214715_mlst.tsv
mlst SRR3215024.fna > SRR3215024_mlst.tsv
mlst SRR3215107.fna > SRR3215107_mlst.tsv
mlst GCF_000009085.1_ASM908v1.fna > GCF_000009085.1_ASM908v1_mlst.tsv

#Rename the first column as requested in the directions and combine into a single file for submission
awk 'BEGIN{FS=OFS="\t"} {sub(/\.fna$/, "", $1)} 1' SRR3214715_mlst.tsv > SRR3214715_mlst_clean.tsv
awk 'BEGIN{FS=OFS="\t"} {sub(/\.fna$/, "", $1)} 1' SRR3215024_mlst.tsv > SRR3215024_mlst_clean.tsv
awk 'BEGIN{FS=OFS="\t"} {sub(/\.fna$/, "", $1)} 1' SRR3215107_mlst.tsv > SRR3215107_mlst_clean.tsv
awk 'BEGIN{FS=OFS="\t"} {sub(/\.fna$/, "", $1)} 1' GCF_000009085.1_ASM908v1_mlst.tsv > GCF_000009085.1_ASM908v1_mlst_clean.tsv
(cat SRR3214715_mlst_clean.tsv && cat SRR3215024_mlst_clean.tsv && cat SRR3215107_mlst_clean.tsv && cat GCF_000009085.1_ASM908v1_mlst_clean.tsv) > mlst.tsv

#Switch working directory, environment to one with checkm, and copy files over for ease (they're small enough)
cd ~/hw3_gtqa/checkm/asm
cp -a ~/hw3_gtqa/fastani/SRR3214715.fna ~/hw3_gtqa/checkm/asm
conda deactivate
conda activate checkm

#Download the database and make sure checkm can find it
cd ~/hw3_gtqa/checkm/db
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz
tar zxvf checkm_data_2015_01_16.tar.gz
echo 'export CHECKM_DATA_PATH=$HOME/hw3_gtqa/checkm/db' >> ~/.bashrc
source ~/.bashrc
echo "${CHECKM_DATA_PATH}"
conda activate checkm

#Run checkm to determine completeness and contamination level
cd ~/hw3_gtqa/checkm
checkm taxon_list | grep Campylo
checkm taxon_set species "Campylobacter jejuni" Cj.markers
checkm analyze Cj.markers ~/hw3_gtqa/checkm/asm analyze_output
checkm qa -f checkm.tax.qa.out -o 1 Cj.markers analyze_output
sed 's/ \+ /\t/g' checkm.tax.qa.out > checkm.tax.qa.out.tsv
cut -f 2- checkm.tax.qa.out.tsv > tmp.tab && mv tmp.tab checkm.tax.qa.out.tsv
sed -i '1d; 3d; $d' checkm.tax.qa.out.tsv

#Rename for submission
cp checkm.tax.qa.out.tsv quality.tsv

#Create a compressed tarball for submission to Canvas
cd ~/hw3_gtqa
tar -czvf assembly_assessment.tar.gz cmds.sh -C ~/hw3_gtqa/fastani fastani.tsv -C ~/hw3_gtqa/mlst mlst.tsv -C ~/hw3_gtqa/checkm quality.tsv