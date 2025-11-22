#!/bin/bash
#SBATCH -n 4
#SBATCH -p main
#SBATCH --mem=100G

# NP Hofford Feb 10 2024

# This is the pipeline script for TE analysis consisting of 6 steps (1) FastQC which analyzes raw sequence reads for quality, (2) Trimmomatic which removes adapter sequences
# (3) Bowtie which removes organellar, bacterial and fungal sequences, (4) Nubeemdedup whoch removes PCR duplicates, (5) interleaving PE sequences for transposome analysis,
# (6) Transposome which analyzes TE content and is called using the multi_transposome.pl script.

# This script is called from the 1_callTEpipe.pl script in the same directory.


# Make overarching directory for sample
mkdir $3
cd $3


##########
# STEP 1 #
##########
# FastQC on raw reads
# Load modules
module load java/11.0.2
# Make FastQC directory
mkdir 1_FastQC
cd 1_FastQC
# FastQC call - $1 and $2 are the R1 and R2 raw reads provided in the initial call
/home/nphofford/TE_pipeline/Installers/FastQC/fastqc $1 $2 -o ./
cd ../


##########
# STEP 2 #
##########
# Adapter trim reads and standardize length
# Make trimmed reads directory
mkdir 2_TrimmedReads
cd 2_TrimmedReads
# Trimmomatic path variables
TRIM_PATH="/mrm/bin/Fast-Plast_JG/Fast-Plast/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
ILLUMINA="/mrm/bin/Fast-Plast_JG/Fast-Plast/bin/adapters/NEB-PE.fa"
# Trimmomatic call - $1 and $2 are raw reads provided in the initial call
java -jar "$TRIM_PATH" PE -threads 4 -trimlog trimlog.txt $1 $2 $3\.trimmed_P1.fq $3\.trimmed_U1.fq $3\.trimmed_P2.fq $3\.trimmed_U2.fq ILLUMINACLIP:"$ILLUMINA":1:30:10 SLIDINGWINDOW:10:20 MINLEN:140
cd ../


##########
# STEP 3 #
##########
# Removal of organellar and contaminated reads
# Make overall bowtie directory
mkdir 3_Bowtie
cd 3_Bowtie
# Load Bowtie2
module load bio/bowtie/2.3

# Make bowtie directory for organellar reads removal
mkdir 3a_Organellar
cd 3a_Organellar
# Bowtie call to remove organellar reads
bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass.fq --un-conc $3\.PE_pass.fq --al $3\.UP_fail.fq --al-conc $3\.PE_fail.fq -p 1 -x /mrm/Sorghum_halepense_Utilities/Organellar/Sorghum_organellar_genomes -1 ../../2_TrimmedReads/$3\.trimmed_P1.fq -2 ../../2_TrimmedReads/$3\.trimmed_P2.fq -S $3.sam
cd ../

# Make bowtie directory for fungal contamination removal
#mkdir 3b_Fungal
#cd 3b_Fungal
# Bowtie call to remove fungal contaminated read
#bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass_fun.fq --un-conc $3\.PE_pass_fun.fq --al $3\.UP_fail_fun.fq --al-conc $3\.PE_fail_fun.fq -p 1 -x /scratch/mrmckain/TE_pipeline/Bowtie_data/Fungal_Test/Fungal_bowtie -1 ../3a_Organellar/$3\.PE_pass.1.fq -2 ../3a_Organellar/$3\.PE_pass.2.fq -S $3\_fun.sam
#cd ../

# EDITED to go directly from the organellar filtering - fungal database not set up
# Make bowtie directory for bacterial contamination removal
mkdir 3c_Bacterial
cd 3c_Bacterial
# Bowtie call to remve bacterial contamination
bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass_clean.fq --un-conc $3\.PE_pass_clean.fq --al $3\.UP_fail_clean.fq --al-conc $3\.PE_fail_clean.fq -p 1 -x /mrm/bin/Bacterial_Genome_DB/NCBI_Bacterial_bowtie -1 ../3a_Organellar/$3\.PE_pass.1.fq -2 ../3a_Organellar/$3\.PE_pass.2.fq -S $3\_bact.sam
cd ../../



##########
# STEP 4 #
##########
# Deduplication of decontaminated reads
# Make deduplication directory
mkdir 4_Deduplication
cd 4_Deduplication
# Nubeamdedup call
# Note that '-r 1' indicates I want the removed duplicates saved to files
/mrm/bin/nubeamdedup-master/Linux/nubeam-dedup -i1 ../3_Bowtie/3c_Bacterial/$3\.PE_pass_clean.1.fq -i2 ../3_Bowtie/3c_Bacterial/$3\.PE_pass_clean.2.fq -r 1
cd ../


##########
# STEP 5 #
##########
# Interleave files for running on Transposome
# Make interleaved directory
mkdir 5_Interleaved
cd 5_Interleaved
# Copy both decontaminated/non-deduplicated and decontaminated/deduplicated fastq files
# We are interested in checking the effect of deduplicated on the transposome run
# Path to decontaminated/deduplicated reads
cp ../4_Deduplication/*uniq.fastq .
# Path to decontaminated/non-deduplicated reads
cp ../3_Bowtie/3c_Bacterial/$3\.PE_pass_clean.*.fq .
# Rename deduplicated files from fastq to fq
mv $3\.PE_pass_clean.1.uniq.fastq $3\.PE_pass_clean.1.uniq.fq
mv $3\.PE_pass_clean.2.uniq.fastq $3\.PE_pass_clean.2.uniq.fq
# Non-deduplicated reads
perl /mrm/bin/fastq2fasta_interleaved.pl $3\.PE_pass_clean.1.fq $3\.PE_pass_clean.2.fq
# Deduplicated reads
perl /mrm/bin/fastq2fasta_interleaved.pl $3\.PE_pass_clean.1.uniq.fq $3\.PE_pass_clean.2.uniq.fq

rm *fq
cd ../



##########
# STEP 6 #
##########
# Run Transposome on prepped reads
# Make Transposome overarching directory
mkdir 6_Transposome
cd 6_Transposome
# Directory for cleaned and deduplicated transposome outut
mkdir 6a_CleanDedup
cd 6a_CleanDedup
# Copy interleaved file to transposome directory
cp ../../5_Interleaved/$3\.PE_pass_clean.1.uniq.fasta .
# Transposome multiple runs call
perl /home/nphofford/TE_pipeline/3_multiTransposome.pl $3\.PE_pass_clean.1.uniq.fasta $3
cd ../

######################
# STEP 6b - Optional #
######################
# Transposome on non-deduplicated reads to check
# Directory for cleaned but NOT deduplicated transposome output
mkdir 6b_CleanNODedup
cd 6b_CleanNODedup
# Copy interleaved file to directory
cp ../../5_Interleaved/$3\.PE_pass_clean.1.fasta .
# Transposome call
perl /home/nphofford/TE_pipeline/3_multiTransposome.pl $3\.PE_pass_clean.1.fasta $3


# DONE
echo "Nate is bad at basketball. LeBron told me.";


