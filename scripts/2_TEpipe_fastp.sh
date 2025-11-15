#!/bin/bash
#SBATCH -n 4
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
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
# Adapter trim reads and standardize length
# Make trimmed reads directory
mkdir 1_TrimmedReads
cd 1_TrimmedReads
# Fastp call - $1 and $2 are raw reads provided in the initial call
# Run using -f (remove first 5 forward reads), -D (deduplicates PCR duplicates), -g (removes polyg tail)
/home/nphofford/fastp/fastp -i $1 -o $3\.R1_filtered.fq -I $2 -O $3\.R2_filtered.fq -h ./$3\_fastp.html -w 4 -f 5 -D -g
cd ../


##########
# STEP 2 #
##########
# Removal of organellar and contaminated reads
# Make overall bowtie directory
mkdir 2_Bowtie
cd 2_Bowtie
# Load Bowtie2
module load bio/bowtie/2.3

# Make bowtie directory for organellar reads removal
mkdir 2a_Organellar
cd 2a_Organellar
# Bowtie call to remove organellar reads
bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass.fq --un-conc $3\.PE_pass.fq --al $3\.UP_fail.fq --al-conc $3\.PE_fail.fq -p 1 -x /mrm/Sorghum_halepense_Utilities/Organellar/Sorghum_organellar_genomes -1 ../../1_TrimmedReads/$3\.R1_filtered.fq -2 ../../1_TrimmedReads/$3\.R2_filtered.fq -S $3.sam
cd ../

# DEPRECTATED
# Make bowtie directory for fungal contamination removal
#mkdir 3b_Fungal
#cd 3b_Fungal
# Bowtie call to remove fungal contaminated read
#bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass_fun.fq --un-conc $3\.PE_pass_fun.fq --al $3\.UP_fail_fun.fq --al-conc $3\.PE_fail_fun.fq -p 1 -x /scratch/mrmckain/TE_pipeline/Bowtie_data/Fungal_Test/Fungal_bowtie -1 ../3a_Organellar/$3\.PE_pass.1.fq -2 ../3a_Organellar/$3\.PE_pass.2.fq -S $3\_fun.sam
#cd ../

# EDITED to go directly from the organellar filtering - fungal database not set up
# Make bowtie directory for bacterial contamination removal
mkdir 2c_Bacterial
cd 2c_Bacterial
# Bowtie call to remve bacterial contamination
bowtie2 --very-sensitive-local --quiet --un $3\.UP_pass_clean.fq --un-conc $3\.PE_pass_clean.fq --al $3\.UP_fail_clean.fq --al-conc $3\.PE_fail_clean.fq -p 1 -x /mrm/bin/Bacterial_Genome_DB/NCBI_Bacterial_bowtie -1 ../2a_Organellar/$3\.PE_pass.1.fq -2 ../2a_Organellar/$3\.PE_pass.2.fq -S $3\_bact.sam
cd ../../



##########
# STEP 3 #
##########
# Interleave files for running on Transposome
# Make interleaved directory
mkdir 3_Interleaved
cd 3_Interleaved
# Path to decontaminated/non-deduplicated reads
cp ../2_Bowtie/2c_Bacterial/$3\.PE_pass_clean.*.fq .
# Non-deduplicated reads
perl ../../scripts/fastq2fasta_interleaved.pl $3\.PE_pass_clean.1.fq $3\.PE_pass_clean.2.fq

rm *fq
cd ../



##########
# STEP 4 #
##########
# Run Transposome on prepped reads
# Make Transposome overarching directory
mkdir 4_Transposome
cd 4_Transposome
# Copy interleaved file to directory
cp ../3_Interleaved/$3\.PE_pass_clean.1.fasta .
# Transposome call
perl ../../scripts/3_multiTransposome.pl $3\.PE_pass_clean.1.fasta $3


# DONE
echo "Nate is bad at basketball. LeBron told me.";


