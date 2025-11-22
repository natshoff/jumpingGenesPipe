# jumpingGenesPipe

A bioinformatics pipeline for estimating transposable element (TE) composition from low-coverage shotgun sequencing data.

## Overview

This pipeline processes Illumina paired-end reads (150 bp) to generate transposable element estimates for *Sorghum halepense*. The workflow includes quality control, decontamination, and TE annotation using Transposome with bootstrapping for robust estimates.

**Study Goal**: Characterizing intra- and inter-population genomic variation through TE composition analysis.

## Pipeline Workflow

The pipeline consists of four main stages:

### Stage 1: Sample Job Submission (`1_callTEpipe.pl`)
- Reads a sample list file
- Locates paired-end FASTQ files for each sample
- Submits batch jobs to SLURM scheduler

### Stage 2: Quality Control & Decontamination (`2_TEpipe_fastp.sh`)
Processes each sample through four steps:

**Step 1 - Quality Filtering & Adapter Trimming (fastp)**
- Removes adapters and low-quality bases
- Trims first 5 bases from forward reads
- Removes PCR duplicates
- Removes polyG tails (common in NovaSeq data)

**Step 2 - Read Decontamination (Bowtie2)**
- **2a**: Removes organellar reads (chloroplast/mitochondria)
- **2c**: Removes bacterial contamination
- Retains only nuclear genomic reads

**Step 3 - Read Interleaving**
- Converts paired-end FASTQ to interleaved FASTA format
- Prepares reads for Transposome analysis

**Step 4 - TE Analysis Setup**
- Initiates Transposome analysis via `3_multiTransposome.pl`

### Stage 3: Bootstrapped TE Estimation (`3_multiTransposome.pl`)
- Generates multiple random subsamples from cleaned reads
- **Default**: 100 iterations of 200,000 read subsamples
- Creates configuration files for each iteration
- Submits individual Transposome jobs

### Stage 4: TE Annotation (`4_transposome.srun`)
- Runs Transposome on each subsample
- Performs BLAST against grass TE database
- Clusters similar sequences
- Annotates repeat families

## Requirements

### Software Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| Perl | 5.x | Pipeline scripts |
| fastp | [version] | Quality control & adapter trimming |
| Bowtie2 | 2.3 | Contamination removal |
| Transposome | 0.12.1 | TE annotation |
| SLURM | - | Job scheduling |

### Reference Databases

1. **Sorghum organellar genomes** (Bowtie2 index)
   - Location: `/mrm/Sorghum_halepense_Utilities/Organellar/Sorghum_organellar_genomes`

2. **NCBI Bacterial genomes** (Bowtie2 index)
   - Location: `/mrm/bin/Bacterial_Genome_DB/NCBI_Bacterial_bowtie`

3. **Grass TE database** (RepBase 21.10)
   - Location: `/mrm/RepBase21.10.fasta/grasrep.ref`
   - Required for TE annotation in `GRStransposome_config.yml`

### System Requirements

- **Memory**: 100GB per sample (main pipeline)
- **CPUs**: 4 cores per sample
- **Storage**: Sufficient scratch space for intermediate files

## Directory Structure

```
jumpingGenesPipe/
├── scripts/
│   ├── 1_callTEpipe.pl              # Job submission script
│   ├── 2_TEpipe_fastp.sh            # Main processing pipeline
│   ├── 3_multiTransposome.pl         # Transposome bootstrapping
│   ├── 4_transposome.srun            # Transposome execution
│   ├── fastq2fasta_interleaved.pl    # Format conversion
│   ├── GRStransposome_config.yml     # Transposome configuration
│   ├── summary.srun                  # Results summarization
│   └── WORKSpull_annotations_transposome_average.pl  # Parse results
├── samples/
│   └── samples.txt                   # List of sample IDs
├── data/
│   ├── raw/                          # Raw FASTQ files
│   └── processed/                    # Processed sample directories
│       └── SAMPLE_NAME/
│           ├── 1_TrimmedReads/
│           ├── 2_Bowtie/
│           │   ├── 2a_Organellar/
│           │   └── 2c_Bacterial/
│           ├── 3_Interleaved/
│           └── 4_Transposome/
└── logs/
    └── slurm-*.out                   # SLURM log files
```

## Usage

### 1. Setup

Clone the repository to your HPC scratch directory:

```bash
cd /scratch/YOUR_USERNAME
git clone git@github.com:natshoff/jumpingGenesPipe.git
cd jumpingGenesPipe
```

### 2. Prepare Sample List

Create a text file with sample identifiers (one per line):

```bash
# samples/samples.txt
JG_CO_011_C
JG_CO_099_A
JG_CO_150_B
```

### 3. Configure File Paths

Edit `scripts/1_callTEpipe.pl` to set your raw data location:

```perl
# Line 23: Update to your sequence directory
my @files = grep { !/\.html$/ } </path/to/your/sequences/*$tarray[0]*R[1-2]*>;
```

### 4. Run Pipeline

Submit jobs for all samples:

```bash
# Run full pipeline (all steps)
perl scripts/1_callTEpipe.pl samples/samples.txt scripts/2_TEpipe_fastp.sh

# Or start from a specific step (e.g., Step 3)
perl scripts/1_callTEpipe.pl samples/samples.txt scripts/2_TEpipe_fastp.sh 3
```

**Available steps:**
- Step 1: Quality Filtering & Adapter Trimming
- Step 2: Contamination Removal  
- Step 3: Read Interleaving
- Step 4: Transposome Analysis

This will submit one SLURM job per sample.

### 5. Monitor Progress

Check job status:

```bash
# View active jobs
squeue -u $USER

# Monitor a specific job log
tail -f logs/slurm-<jobid>.out
```

### 6. Restart Failed Jobs (Optional)

If a job fails at a specific step, you can restart from that step without reprocessing earlier steps:

```bash
# Create a sample list with only the failed samples
echo "JG_CO_011_C" > samples/failed_samples.txt

# Restart from Step 3 (skipping Steps 1-2)
perl scripts/1_callTEpipe.pl samples/failed_samples.txt scripts/2_TEpipe_fastp.sh 3
```

The pipeline will skip completed steps and use existing intermediate files.

### 7. Summarize Results

**[WIP]** - Results summarization documentation coming soon.

## Configuration

### Adjusting Transposome Sampling Parameters

In `scripts/3_multiTransposome.pl`, you can modify:

```perl
# Line 14: Number of bootstrap iterations (default: 100)
for (my $i=0; $i<=99; $i++){

# Line 19 & 22: Number of reads per subsample (default: 200,000)
$start_rand = int(rand($total_pairs-200000));
`tail -n $start_rand $ARGV[0] | head -n 200000 > ...`;
```

### Modifying Transposome Parameters

Edit `scripts/GRStransposome_config.yml`:

```yaml
clustering_options:
  - percent_identity:   90      # Clustering threshold
  - fraction_coverage:  0.55    # Minimum alignment coverage

annotation_options:
  - cluster_size:       100     # Minimum cluster size for annotation
```

See the [Transposome documentation](https://github.com/sestaton/Transposome/wiki) for more details.


## Output Files

**[WIP]** - Detailed output documentation coming soon.

Key output locations:
- Cleaned reads: `data/processed/SAMPLE/2_Bowtie/2c_Bacterial/`
- Transposome results: `data/processed/SAMPLE/4_Transposome/*_transposome_*_out/`
- TE annotations: `*_report_annotations_summary.tsv`

## Troubleshooting

### Log Files

- **Pipeline logs**: `logs/slurm-<jobid>.out`
- **Transposome logs**: `data/processed/SAMPLE/4_Transposome/*_log.txt`

## Citation

If you use this pipeline, please cite:

**Transposome:**
- Staton SE, Burke JM (2015). Transposome: A toolkit for annotation of transposable element families from unassembled sequence reads. *Bioinformatics* 31(11): 1827-1829.

**fastp:**
- Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17): i884-i890.

**Bowtie2:**
- Langmead B, Salzberg SL (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9(4): 357-359.

## Author

NP Hofford  
University of Alabama

## License

[Add license information]

## Acknowledgments

- McKain Lab for computational resources
- RepBase for TE reference database

